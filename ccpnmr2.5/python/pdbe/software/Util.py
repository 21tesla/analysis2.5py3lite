
"""
======================COPYRIGHT/LICENSE START==========================

Util.py: Useful functions for software in this directory

Copyright (C) 2008 Wim Vranken (European Bioinformatics Institute)

=======================================================================

Can only be distributed and used within the CCPNGRID project, no other
use or distribution allowed.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
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

#
# General functions
#


def getResMapping(resonances, ignoreSerials = None, ignoreChains = None, ignoreResidues = None, onlyElements = None):

  """
  
  Function that maps resonances to residues and atom names
  
  """

  if not ignoreSerials:
    ignoreSerials = []
    
  if not ignoreChains:
    ignoreChains = []

  if not ignoreResidues:
    ignoreResidues = {}

  if not onlyElements:
    onlyElements = []

  resMapping = {}
  assignedResidues = []
  excludedResonances = []
  
  #
  # First make quick link for resonance -> atom
  #
  
  for resonance in resonances:

    if resonance.serial in ignoreSerials:
      continue

    #
    # Only a link from the resonance to an atom if there is a resonanceSet...
    #
    
    if resonance.resonanceSet:
    
      atomSets = list(resonance.resonanceSet.atomSets)
      residue = atomSets[0].findFirstAtom().residue
      chainCode = residue.chain.code
      
      if chainCode in ignoreChains:
        excludedResonances.append(resonance)
        continue
            
      if ignoreResidues.has_key(chainCode) and residue.seqId in ignoreResidues[chainCode]:
        excludedResonances.append(resonance)
        continue
        
      #
      # Go over the atomSets...
      #
      
      atomNameList = []

      for atomSet in atomSets:
        
        refAtom = atomSet.findFirstAtom()
        curResidue = refAtom.residue
        
        #
        # Check if all is OK (no mapping to different residues)
        #

        if curResidue != residue:
          print "  ERROR two residues to same resonance!"
          atomNameList = []
          break
          
        #
        # Include only the element types that are in the onlyElements list
        # (if empty, will include all)
        #
        
        if onlyElements and refAtom.chemAtom.elementSymbol not in onlyElements:
          atomNameList = []
          break

        atomNameList.append(atomSet.name)
      
      if atomNameList:
      
        atomNameList.sort()
        atomNameTuple = tuple(atomNameList)

        resMapping[resonance] = [residue,atomNameTuple]

        if not residue in assignedResidues:
          assignedResidues.append(residue)
      
      else:
        excludedResonances.append(resonance)

    
    else:
    
      #
      # If the resonance is not linked and there is an onlyElements list provided,
      # check whether it has a valid name, and whether its atom name falls within
      # the onlyElements list (if available).
      #
    
      if onlyElements:
        
        resNames = getApplResNames('nmrStar',[resonance])
        
        if resNames:
        
          resName = resNames.keys()[0] # Should be sufficient at this stage

          (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName, verbose = 0)

          if chainCode == seqCode == atomName == None:
            excludedResonances.append(resonance)
            print "  Removed invalid resonance."

          elif atomName:

            # This is not 100% foolproof but should be helpful

            elementPresent = False

            for elementSymbol in onlyElements:
              if elementSymbol in atomName:
                elementPresent = True
                break
                
            if not elementPresent:
              excludedResonances.append(resonance)
              print "  Removed resonance %s - not connected and bad element." % resName

  return (resMapping,assignedResidues,excludedResonances)

#
# General HTML setup
#

import os, time

class HtmlPage:

  """
  Definition of HTML page. Origin from pdbe/web/Util.py and pdbe/coco/Util.py HtmlPage definitions
  Can define own stylesheet.
  """

  def __init__(self,fileName, styleSheet = "", htmlBaseName = None):
  
    self.fileName = fileName    
    self.htmlText = ""
    
    self.styleSheet = styleSheet

    self.initialize()

  def setupHtml(self,title,colspan = 2):

    self.write("""
  <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
  <html>
  <head>
  <title>%s</title>
  <link rel="stylesheet" type="text/css" href="%s" title="Style">
  <script src="http://www.ebi.ac.uk/pdbe/search/help.js" type="text/javascript"></script>
  </head>

  <body>   
  <table width="100%%">
    <tr><td height="50" valign="middle" colspan=%d><h1>%s</h1></td></tr>
    """ % (title,self.styleSheet,colspan,title))
    
  def setupSideMenu(self,subMenuItems):
    
    pass
      
  def createFullSideMenu(self,fullSideMenuList, indent = "    ", activeItem = None, currentItem = tuple()):
    
    pass
                 
  def mainPageHtml(self,titleString):
    
    self.writeTableHeader(titleString, colspan = 99)

  def addMainTableColumn(self,text,addText = ''):
  
    self.write('        <td%s>%s</td>\n' % (addText,text))

  def writeTableHeader(self,tableHeader,colspan = 2):
  
    self.addMainTableRow()
    self.write('<th colspan="%d">%s</th>' % (colspan,tableHeader))
    self.closeMainTableRow()

  def addSubHeading(self,subHeadingText, colspan = 2):
  
    self.addMainTableRow()
    self.addMainTableColumn(subHeadingText, addText = ' colspan="%d" class = "subheading"' % (colspan))
    self.closeMainTableRow()

  def finishHtml(self,creationString = None, timeString = None, developerText = None, hasMainTable = False):
  
    if not timeString:
      timeString = time.strftime("%a %b %d %H:%M:%S GMT %Y", time.gmtime(time.time()))

    if not creationString:
      creationString = timeString
      
    if hasMainTable:
      tableText = "</table>"
    else:
      tableText = ""
    
    self.write("""
<div class="banner">
<table width="100%%">
%s

	 <td align="right">
  <!-- Created: %s -->
  Last modified: %s
	  </td>
  </tr>
</table>
</div>
%s

</body>
</html>
    """ % (self.bottomBar(developerText = developerText),creationString,timeString,tableText))
    
    self.writeHtmlFile()

  def bottomBar(self,developerText = None):
  
    if not developerText:
      developerText = '<a href="mailto:wim@ebi.ac.uk">Wim Vranken</a>'
  
    bottomBar = """
      <tr valign="top">
      <td>Primary developers: %s</td>
    """ % developerText
    
    return bottomBar
    
  def initialize(self):
      
    self.setDirs()
    
  def setDirs(self):
  
    (dirName,fileName) = os.path.split(self.fileName)
    fullDirName = os.path.join(os.getcwd(),dirName)
    
    self.createDir(fullDirName)
      
  def createDir(self,fullDirName):
  
    if not os.path.exists(fullDirName):
      (baseDirName,dirName) = os.path.split(fullDirName)
      if not os.path.exists(baseDirName):
        self.createDir(baseDirName)
      os.mkdir(fullDirName)
 
  def setInfoValueRow(self,info,value,infoAddText = '', valueAddText = ' class = "centered"'):
  
    info = info.replace(" ","&nbsp;")
    info = info.replace("<a&nbsp;","<a ")
    value = value.replace(" ","&nbsp;")
  
    self.addMainTableRow()
    self.addMainTableColumn(info, addText = ' class="leftsubheading" %s' % infoAddText) 
    self.addMainTableColumn(value, addText = '%s' % valueAddText)
    self.closeMainTableRow()
    
  def setEmptyRow(self,colspan = 2):
  
    self.addMainTableRow()
    self.addMainTableColumn('&nbsp;', addText = ' colspan="%d"' % colspan)
    self.closeMainTableRow()
  
  def fillColumns(self,fromCol,toCol):
  
    for colNum in range(fromCol,toCol,1):
      self.addMainTableColumn('&nbsp;')

  def write(self,htmlText):
  
    self.htmlText += htmlText
    
  def writeHtmlFile(self):
  
    outFile = open(self.fileName,'w')
    outFile.write(self.htmlText)
    outFile.close()

  def addMainTableRow(self):
  
    self.write("      <tr>\n")
  
  def addMainTableColumn(self,text,addText = ''):
  
    self.write('        <td valign="middle"%s>%s</td>\n' % (addText,text))

  def closeMainTableRow(self):
  
    self.write("      </tr>\n")

  def closeMainTable(self):
  
    self.write('	     </table>\n')
    
    
#
# Handling of contactShortDict information to get occurrences for interatom contacts
#

from ccp.general.Util import getSortContactDist

class ContactOccurrenceHandler:

  contactShortDict = {}
  maxRangeCheckDict = {}
  residueSecStrucDict = {}

  distanceClassesDict = {}
  distanceClassesInfoDict = {}

  secStrucKeysDict = {}
  secStrucKeysInfoDict = {}
  
  resonancesStore = None  # Can be either NmrConstraintStore or NmrProject
  
  contactInfoImportPath = "pdbe.software.data"
  
  class ContactOccurrenceError(StandardError):
    
    def __init__(self, value):
      self.value = value
    def __str__(self):
      return repr(self.value)
   
  def loadCsdModule(self,moduleString):
  
    try:
      fullModule = '%s.%s' % (self.contactInfoImportPath,moduleString)
      contactInfoModule = __import__(fullModule,{},{},moduleString)
      contactInfo = getattr(contactInfoModule,'contactShortDict')
      referenceInfo = getattr(contactInfoModule,'referenceInfo')
      print 'Imported %s' % fullModule
      
    except:
      contactInfo = None
      referenceInfo = None
      print 'Could not import dictionary %s - no information available' % moduleString
      
    return (contactInfo,referenceInfo)

  def importContactShortDicts(self):
  
    chainCodes = self.chainShortInfo.keys()
    chainCodes.sort()
    
    chainCodeCombs = []
    
    for chainCode in chainCodes:
      for otherChainCode in chainCodes:
        chainCodeCombs.append((chainCode,otherChainCode))
   
    csdModuleLoaded = {}
    for (chainCode,otherChainCode) in chainCodeCombs:
    
      csiInfo = []

      for tmpChainCode in (chainCode,otherChainCode):
        csiInfo.append(self.chainShortInfo[tmpChainCode])
      csiInfo.sort()

      if csiInfo[0] == csiInfo[1]:
        # Same, easy
        if csiInfo[0][1] and self.useProteinFamilies:
          csdKey = "%s_%s" % (csiInfo[0][0],csiInfo[0][1])
        else:
          csdKey = "%s_all" % (csiInfo[0][0])

      elif csiInfo[0][0] == csiInfo[1][0]:
        # Protein but different type
        protFams = [csiInfo[0][1],csiInfo[1][1]]
        protFams.sort()

        # TODO should probably have a separate class for intermolecular... anyway just stick to alpha_beta if possible
        if not protFams.count(None) and self.useProteinFamilies:
          if protFams[0] == protFams[1]:
            protFamKey = protFams[0]
          elif csiInfo[0][0] == 'protein':
            protFamKey = 'alpha_beta'
          else:
            protFamKey = "_".join(protFams)

          csdKey = "%s_%s" % (csiInfo[0][0],protFamKey)

        else:
          if protFams.count(None) == 1 and self.useProteinFamilies:
            protFam = protFams[not protFams.index(None)]
            csdKey = "%s_%s" % (csiInfo[0][0],protFam)
          else:
            csdKey = "%s_all" % (csiInfo[0][0])

      else:
        csdKeyEls = []
        for (molType,protFam) in csiInfo:
          if not protFam or not self.useProteinFamilies:
            csdKeyEls.append("%s_all" % (molType))
          else:
            csdKeyEls.append("%s_%s" % (molType,protFam))

        csdKey = "_".join(csdKeyEls)
        
      
      if csdKey not in csdModuleLoaded.keys():      
        csdModuleLoaded[csdKey] = self.loadCsdModule("contactInfo_%s" % csdKey)

      self.contactShortDict[(chainCode,otherChainCode)] = csdModuleLoaded[csdKey][0]
      self.maxRangeCheckDict[(chainCode,otherChainCode)] = csdModuleLoaded[csdKey][1]['maxRangeCheck']
      
      self.secStrucKeysDict[(chainCode,otherChainCode)] = csdModuleLoaded[csdKey][1]['secStrucKeys']
      self.secStrucKeysInfoDict[(chainCode,otherChainCode)] = {'number': len(self.secStrucKeysDict[(chainCode,otherChainCode)])}
      
      self.distanceClassesDict[(chainCode,otherChainCode)] = csdModuleLoaded[csdKey][1]['distanceClasses']
      
      # Precalculate distance step, lower and upper cutoff 
      distStep = self.distanceClassesDict[(chainCode,otherChainCode)][1] - self.distanceClassesDict[(chainCode,otherChainCode)][0]
      self.distanceClassesInfoDict[(chainCode,otherChainCode)] = {'step': distStep,
                                                                  'number': len(self.distanceClassesDict[(chainCode,otherChainCode)]),
                                                                  'medianIndex': int(len(self.distanceClassesDict[(chainCode,otherChainCode)]) / 2.0),
                                                                  'lowerCutoff': self.distanceClassesDict[(chainCode,otherChainCode)][0] - distStep,
                                                                  'upperCutoff': self.distanceClassesDict[(chainCode,otherChainCode)][-2]}  # -1 is the 'long' class!
        
  def setInvConfidence(self,confidence):
  
    if 0.0 <= confidence <= 1.0:
      self.invConfidence = 1.0 - confidence
    else:
      print "  Confidence has to be between 0 and 1 - resetting to default %.2f" % self.defaultConfidence
      self.invConfidence = 1.0 - self.defaultConfidence
          
  def createAtomSetInfo(self,project):
    
    #
    # Get info for atomsets, in case single atoms have to be 'expanded' to multiple ones.
    # Should only be the case for prochirals (????)
    #
    
    self.atomSetInfo = {}
    self.prochiralAtomInfo = {}
    
    for molType in self.molTypes:
      self.atomSetInfo[molType] = {}
      for ccpCode in self.ccpCodes[molType]:
        self.atomSetInfo[molType][ccpCode] = {}
        
        chemComp = project.findFirstChemComp(molType = molType, ccpCode = ccpCode)
        
        for chemAtomSet in chemComp.findAllChemAtomSets(isProchiral = True):

          if chemAtomSet.chemAtoms:
            chemAtomOrSets = chemAtomSet.chemAtoms
          else:
            chemAtomOrSets = chemAtomSet.chemAtomSets
          
          # Also track prochiral atoms - in case need alternative
          if chemAtomSet.isProchiral:
            if not self.prochiralAtomInfo.has_key(molType):
              self.prochiralAtomInfo[molType] = {}
            if not self.prochiralAtomInfo[molType].has_key(ccpCode):
              self.prochiralAtomInfo[molType][ccpCode] = {}
              
          for chemAtomOrSet in chemAtomOrSets:
            chemAtomOrSetNames = [chemAtomOrSet.name for chemAtomOrSet in chemAtomOrSets]
            chemAtomOrSetNames.sort()
            chemAtomOrSetNamesTuple = tuple(chemAtomOrSetNames)
            for chemAtomOrSetName in chemAtomOrSetNames:
              self.atomSetInfo[molType][ccpCode][chemAtomOrSetName] = chemAtomOrSetNamesTuple
              
              if chemAtomSet.isProchiral:
                self.prochiralAtomInfo[molType][ccpCode][chemAtomOrSetName] = chemAtomOrSetNames[not chemAtomOrSetNames.index(chemAtomOrSetName)]
              
    
    #
    # Custom info - can I easily set this from ref data?
    #
    
    molType = 'protein'
    if molType in self.molTypes:
      for ccpCode in ('Phe','Tyr'):
        if ccpCode in self.atomSetInfo[molType]:
          self.atomSetInfo[molType][ccpCode]["CZ"] = ('HE1','HE2')
          self.atomSetInfo[molType][ccpCode]["HE*"] = ('HE1','HE2')
          self.atomSetInfo[molType][ccpCode]["HE1"] = ('HE1','HE2')
          self.atomSetInfo[molType][ccpCode]["HE2"] = ('HE1','HE2')
          self.atomSetInfo[molType][ccpCode]["CG"] = ('HD1','HD2')
          self.atomSetInfo[molType][ccpCode]["HD*"] = ('HD1','HD2')
          self.atomSetInfo[molType][ccpCode]["HD1"] = ('HD1','HD2')
          self.atomSetInfo[molType][ccpCode]["HD2"] = ('HD1','HD2')
      
  def setupContactOccurrenceInfo(self,resonancesStore, hasDistanceInfo=True, residueSecStrucDict=None, useProteinFamilies=False):  
    
    #
    # This is to doublecheck that top-level values don't get mixed with lower level ones..
    #

    self.hasDistanceInfo = hasDistanceInfo
    self.useProteinFamilies = useProteinFamilies

    #
    # Get secondary structure information - this can be reset this way
    #
    
    if residueSecStrucDict:

      self.residueSecStrucDict = residueSecStrucDict
    
    #
    # Set other information if new resonance store
    #

    if self.resonancesStore != resonancesStore:
    
      self.resonancesStore = resonancesStore
      
      #
      # Get the resonances depending on the store type
      #
      
      if self.resonancesStore.className == 'NmrProject':
        resonances = self.resonancesStore.resonances
      else:
        resonances = self.resonancesStore.fixedResonances
      
      #
      # Get resonance to atom mapping info (also on residues, ...)
      #

      (self.resMapping,self.assignedResidues,self.excludedResonances) = getResMapping(resonances)
      
      #
      # Set secondary structure info (if necessary)
      #
  
      self.setResidueSecStrucDict(residueSecStrucDict = self.residueSecStrucDict)
       
      #
      # Determine molecular types and protein family, if possible, for importing dictionary
      #
      
      self.setMolTypesResiduesAndProteinFamily()
                   
      #
      # Automatically import contact occurrence information
      # 

      self.importContactShortDicts()

      #
      # Also create reference info for atomSets
      # 

      self.createAtomSetInfo(resonancesStore.root)
  
  def setupContactOccurrenceInfo_nonCcpnBased(self,molTypes,confidence, hasDistanceInfo=True, useProteinFamilies=False):

    from ccp.general.Constants import standardResidueCcpCodes
    from memops.api import Implementation

    #
    # This is to doublecheck that top-level values don't get mixed with lower level ones..
    #

    self.hasDistanceInfo = hasDistanceInfo
    self.useProteinFamilies = useProteinFamilies

    #
    # Warning: this only sets info for standard residues. Should be good enough,
    # use CCPN if not.
    #

    self.molTypes = molTypes
    self.ccpCodes = {}

    for molType in self.molTypes:
      self.ccpCodes[molType] = standardResidueCcpCodes[molType]

    self.importContactShortDicts()
    self.setInvConfidence(confidence)
    
    project = Implementation.MemopsRoot(name = 'test')
    self.createAtomSetInfo(project)
    
  def setResidueSecStrucDict(self,residueSecStrucDict = None):

    #
    # Get secondary structure information
    #
    
    if residueSecStrucDict:
    
      # Assuming this information is OK!
      
      self.residueSecStrucDict = residueSecStrucDict
      
    else:
    
      from ccp.general.Util import getResidueSsCode
    
      hasSsInfo = False
    
      for residue in self.assignedResidues:
      
        ssCode = getResidueSsCode(residue)
        
        if ssCode:
          hasSsInfo = True
        else:
          ssCode = "C"

        self.residueSecStrucDict[(residue.chain.code,residue.seqId)] = ssCode
        
      #
      # Set to empty if nothing there...
      #
      
      if not hasSsInfo:
        self.residueSecStrucDict = {}

  def setMolTypesResiduesAndProteinFamily(self):
  
    """
    Set the molecular types and ccpCodes per molType, also protein family classes, per chain and per inter-chain (not yet relevant!)
    """
    
    self.molTypes = []
    self.ccpCodes = {}
    self.residueInfo = {None: (None,(None,None),None)}
    
    chainInfo = {}

    # Molecular information
    for residue in self.assignedResidues:

      molType = residue.molType
      ccpCode = residue.ccpCode
      
      chainCode = residue.chain.code
      
      if chainCode not in chainInfo.keys():
        chainInfo[chainCode] = {'molTypes': [], 'secStrucs':[]}    

      chainInfo[chainCode]['molTypes'].append(molType)

      # WARNING: this was changed to make code faster (and contactOccurrence useable from outside CCPN)
      self.residueInfo[residue] = (molType,(residue.chain.code,residue.seqId),ccpCode)

      if not molType in self.molTypes:
        self.molTypes.append(molType)
        self.ccpCodes[molType] = []

      if not ccpCode in self.ccpCodes[molType]:
        self.ccpCodes[molType].append(ccpCode)

    # Now secondary structure
    if self.residueSecStrucDict:
    
      for (chainCode,seqId) in self.residueSecStrucDict.keys():
        
        # No point doing this if there's no assignments anyway...
        if chainCode in chainInfo.keys():
          
          chainInfo[chainCode]['secStrucs'].append(self.residueSecStrucDict[(chainCode,seqId)])
          
    #
    # Cleanup chaincode info...
    #
    
    self.chainShortInfo = {}
    
    
    print chainInfo
    print len(self.assignedResidues)
    
    for chainCode in chainInfo.keys():
      
      self.chainShortInfo[chainCode] = [None,'alpha_beta']
      
      #
      # Set the chain molecular type information
      #
      
      molTypes = chainInfo[chainCode]['molTypes']
      molTypes.sort()
      
      refMolType = molTypes[int(len(molTypes)/2.0)]
      
      if molTypes.count(refMolType) == len(molTypes) or (refMolType != 'other' and molTypes.count(refMolType) + molTypes.count('other') == len(molTypes)):
        self.chainShortInfo[chainCode][0] = refMolType
      else:
        uniqueMolTypes = []
        for molType in molTypes:
          if molType != 'other' and molType not in uniqueMolTypes:
            uniqueMolTypes.append(molType)
        uniqueMolTypes.sort()
        self.chainShortInfo[chainCode][0] = '_'.join(uniqueMolTypes)

      #
      # For proteins, set the family
      #
      
      if self.chainShortInfo[chainCode][0] == 'protein' and chainInfo[chainCode]['secStrucs']:
      
        ssCodes = chainInfo[chainCode]['secStrucs']
        
        protFam = 'alpha_beta'

        total = len(ssCodes)
        fractions = {}
        for (ssType,ssCode) in (('alpha','H'),('coil','C'),('beta','E')):
          fractions[ssType] = ssCodes.count(ssCode) * 1.0 / total
        
        print fractions
        if fractions['alpha'] >= 0.05 and fractions['beta'] >= 0.05:
          protFam = 'alpha_beta'
        elif fractions['beta'] > 0.1:
          protFam = 'beta'
        elif fractions['alpha'] > 0.1:
          protFam = 'alpha'
        elif fractions['coil'] > 0.5:
          protFam = 'coil'
        
        self.chainShortInfo[chainCode][1] = protFam
        
  def convertSsInfoToResidueSecStrucDict(self,molSystem,ssInfo,setMissingToCoil = False):
  
    """
    Function to convert ssInfoDict (with string keys for chains and residues) into
    CCPN object based residue dictionary
    """
  
    residueSecStrucDict = {}

    for chainCode in ssInfo.keys():

      chain = molSystem.findFirstChain(code = chainCode)
      
      if not chain:
        print "  Chain %s missing in secondary structure info conversion, ignored." % chainCode
        continue

      for residueKey in ssInfo[chainCode].keys():

        (seqCode,insertionCode) = residueKey

        residue = chain.findFirstResidue(seqCode = seqCode, seqInsertCode = insertionCode)

        if not residue:
          print "No mapping", chainCode, seqCode

        residueSecStrucDict[(chain.code,residue.seqId)] = ssInfo[chainCode][residueKey]
      
      #
      # Set missing info to coil if required
      #
      
      if setMissingToCoil:
      
        for residue in chain.sortedResidues():
          residueKey = (chain.code,residue.seqId)
          if not residueSecStrucDict.has_key(residueKey):
            residueSecStrucDict[residueKey] = 'C'
      
    return residueSecStrucDict
    
     
  def getContactOccurrence(self,resInfo,debugMode=False, distance=None, ssCodes=None, contactOccurrenceDefault=None):

    """
    Function to get contact occurrence out of a dictionary with finalised results from constraint analysis.
    """
    #
    # Initialise variables
    #

    (molType1,residueKey1,ccpCode1,atomNameTuple1) = resInfo[0]
    (molType2,residueKey2,ccpCode2,atomNameTuple2) = resInfo[1]
    
    (chainCode1,seqId1) = residueKey1
    (chainCode2,seqId2) = residueKey2
    
    chainCodes = [chainCode1,chainCode2]
    chainCodes.sort()
    chainCodesTuple = tuple(chainCodes)

    molTypes = (molType1,molType2)
    
    #
    # Determine the secondary structure codes information
    #
    
    ssCodeTuple = self.getSsCodes(distance,ssCodes,resInfo)
    
    #
    # Order the incoming information
    #
    
    (sortType,contactDist) = getSortContactDist(residueKey1,residueKey2,self.maxRangeCheckDict[chainCodesTuple])
          
    # Don't bother if no info available      
    if not self.contactShortDict.has_key(chainCodesTuple) or not self.contactShortDict[chainCodesTuple]:
      return (contactOccurrenceDefault,None,None)
    else:
      contactShortDict = self.contactShortDict[chainCodesTuple]
      
    # Order info for residues 
    if (sortType == 'seqId' and seqId1 <= seqId2) or \
       (sortType == 'ccpCode' and (ccpCode1 < ccpCode2 or (ccpCode1 == ccpCode2 and atomNameTuple1 < atomNameTuple2))):
      residueNames = (ccpCode1,ccpCode2)
      atomNamesList = [atomNameTuple1,atomNameTuple2]

    else:
      residueNames = (ccpCode2,ccpCode1)
      atomNamesList = [atomNameTuple2,atomNameTuple1]
      ssCodeTuple = (ssCodeTuple[1],ssCodeTuple[0])  # This was missing previously!! Problem!

    # Also sort atom names if same residue (in same chain!)
    if sortType == 'seqId' and seqId1 == seqId2 and chainCode1 == chainCode2:
      atomNamesList.sort()

    #
    # Set up some information for comparing to shortContactDict
    #
    
    (distanceClassIndex,ssCodeIndex) = self.getDistanceClassAndSsCodeIndex(distance,ssCodeTuple,chainCodesTuple)
    
    if debugMode:
      print "RK:",residueNames
      print "CD:",contactDist
      print "AL0:",atomNamesList[0]
      print "AL1:",atomNamesList[1]
      print "ST:", sortType
      print "SS:", ssCodes, ssCodeIndex
      print "DIST:", distance, distanceClassIndex
      print

    #
    # Some more initialisation, start to get the occurrence out
    #

    contactOccurrence = contactOccurrenceDefault
    averageDist = None

    # TODO: if going to use sorted list, would be faster to do this as a list of lists...

    if contactShortDict.has_key(residueNames):
      #if debugMode:
      #  print "  Found residueNames"
      if contactShortDict[residueNames].has_key(contactDist):
        #if debugMode:
        #   print "  Found contactDist"
        
        atomNameKey1 = self.findAtomNameKey(contactShortDict[residueNames][contactDist],molTypes[0],residueNames[0],atomNamesList[0])

        #
        # Continue if match found for first atom name tuple
        #

        if atomNameKey1:
          if debugMode:
            print "SKEY1",atomNameKey1, atomNamesList[0]

          self.residueCcpCode = residueNames[1]
          atomNameKey2 = self.findAtomNameKey(contactShortDict[residueNames][contactDist][atomNameKey1],molTypes[1],residueNames[1],atomNamesList[1])

          if atomNameKey2:
            if debugMode:
              print "SKEY2",atomNameKey2, atomNamesList[1]
              print contactShortDict[residueNames][contactDist][atomNameKey1][atomNameKey2]
              print self.secStrucKeysInfoDict[chainCodesTuple]['number']
            
            #
            # If no distance given or intermolecular, use the overall occurrence - secondary structure ignored.
            # In case of intermolecular, not enough data to get distance dependence
            #
            
            if contactDist == -2 or distanceClassIndex == None:
              contactOccurrence = contactShortDict[residueNames][contactDist][atomNameKey1][atomNameKey2][0] 
              if self.hasDistanceInfo:
                # Normalise in this case!! 'downgrade' the overall occurrence so that matches (more or less) the distance class values
                contactOccurrence = contactOccurrence / self.distanceClassesInfoDict[chainCodesTuple]['number']
            
            else:
              distanceOccurrenceInfo = contactShortDict[residueNames][contactDist][atomNameKey1][atomNameKey2][1][distanceClassIndex]
              
              if distanceOccurrenceInfo[0] == None:
                # No data
                contactOccurrence = 0.0
              
              else:
                # Data available, but is there anything on the ss level?
                if ssCodeIndex == None:
                  # Normalise - don't want this number to be too high!
                  contactOccurrence = distanceOccurrenceInfo[0] / self.secStrucKeysInfoDict[chainCodesTuple]['number']

                else:
                  contactOccurrence = distanceOccurrenceInfo[1][ssCodeIndex]
                  
                  if contactOccurrence == None:
                    # No data, again use the default
                    contactOccurrence = distanceOccurrenceInfo[0] / self.secStrucKeysInfoDict[chainCodesTuple]['number']

    if debugMode and contactOccurrence:
      print residueNames, contactDist, atomNameKey1, atomNameKey2, contactOccurrence

    return (contactOccurrence,atomNamesList,averageDist)
  
  def getSsCodes(self,distance,ssCodes,resInfo):
  
    if not ssCodes:
    
      ssCodes = [None,None]
      
      if self.residueSecStrucDict:
        ssCodes = []
        for residueKey in (resInfo[0][1],resInfo[1][1]):
          ssCode = None
          if self.residueSecStrucDict.has_key(residueKey):
            ssCode = self.residueSecStrucDict[residueKey]
          # Check whether right info is being passed in!!
          elif type(residueKey[0]) != type(self.residueSecStrucDict.keys()[0][0]):
            raise self.ContactOccurrenceError("Problem with secondary structure information, am getting %s but expecting %s!" % (type(self.residueSecStrucDict.keys()[0][0]),type(residueKey[0])))
          ssCodes.append(ssCode)
     
    ssCodeTuple = tuple(ssCodes)
      
    return ssCodeTuple
  
  def getDistanceClassAndSsCodeIndex(self,distance,ssCodeTuple,chainCodesTuple):
    
    distanceClassIndex = ssCodeIndex = None
    
    if distance != None:
      
      if distance >= self.distanceClassesInfoDict[chainCodesTuple]['upperCutoff']:
        distanceClassIndex = -1
      elif distance < self.distanceClassesInfoDict[chainCodesTuple]['lowerCutoff']:
        distanceClassIndex = 0
      else:
        distanceClassIndex = int((distance - self.distanceClassesInfoDict[chainCodesTuple]['lowerCutoff']) / self.distanceClassesInfoDict[chainCodesTuple]['step'])

      # This should always be true, will look for (None,None) ss code set, which is present...
      if ssCodeTuple in self.secStrucKeysDict[chainCodesTuple]:
        ssCodeIndex = self.secStrucKeysDict[chainCodesTuple].index(ssCodeTuple)

    return (distanceClassIndex,ssCodeIndex)

  def findAtomNameKey(self,contactShortDictInfo,molType,ccpCode,atomNameTuple):

    atomNameKey = None

    if contactShortDictInfo.has_key(atomNameTuple):
      atomNameKey = atomNameTuple
    else:
      atomSetTuple = self.getAtomSetTuple(molType,ccpCode,atomNameTuple)

      if atomSetTuple and contactShortDictInfo.has_key(atomSetTuple):
        atomNameKey = atomSetTuple

    return atomNameKey

  def getAtomSetTuple(self,molType,ccpCode,atomNameTuple):

    """
    Gets atom set name from an atom name tuple, used for contact occurrence
    """

    atomSetTuple = None

    if len(atomNameTuple) == 2:
      if atomNameTuple[0][-1] != '*':
        atomSetTuple = (atomNameTuple[0][:-1] + '*',)
      else:
        atomSetTuple = (atomNameTuple[0][:-2] + '*',)

    elif self.atomSetInfo.has_key(molType) and self.atomSetInfo[molType].has_key(ccpCode):
      
      # Should only be one value
      atomName = atomNameTuple[0]

      if self.atomSetInfo[molType][ccpCode].has_key(atomName):

        atomSetTuple = self.atomSetInfo[molType][ccpCode][atomName]

    return atomSetTuple

#
# Class with functions to set up dicts and handle resonance/atom/coordinate stuff
#

class ResonanceCoordinateHandler:

  def setAssignedAtomsAndResidues(self,fixedResonances,chainCodeFilter=None):
  
    from ccp.general.Util import getResAtomObjectMapping
  
    self.assignedResonances = []

    for fres in fixedResonances:
      if fres.resonanceSet:
        self.assignedResonances.append(fres)

    self.resObjectMapping = getResAtomObjectMapping(self.assignedResonances,chainCodeFilter=chainCodeFilter)

    #
    # Track all assigned atoms - no point doing these (in principle)
    #

    self.assignedAtoms = []
    self.assignedResidues = []

    for fres in self.assignedResonances[:]:
      
      if fres not in self.resObjectMapping.keys():
        self.assignedResonances.pop(self.assignedResonances.index(fres))
        continue
      
      (residue,atomList) = self.resObjectMapping[fres]
      for atom in atomList:
        self.assignedAtoms.append(atom)

      if not residue in self.assignedResidues:
        self.assignedResidues.append(residue)
  
  def createCoordAtomInfoDict(self):
  
    self.coordAtomInfo = {}

    for coordChain in self.structureList[0].structureEnsemble.coordChains:
      for coordResidue in coordChain.residues:
        residue = coordResidue.residue
        
        if not residue:
          continue
        
        for atom in residue.atoms:

          if not atom in self.assignedAtoms:
            continue

          if not self.coordAtomInfo.has_key(atom):
            self.coordAtomInfo[atom] = []

          coordAtom = coordResidue.findFirstAtom(atom = atom)
          
          coordsDict = {}
          
          if coordAtom:
            for coord in coordAtom.coords:
              coordsDict[coord.model] = coord
          
          for model in self.structureList:
            coord = None
            if coordsDict.has_key(model):
              coord = coordsDict[model]
            self.coordAtomInfo[atom].append(coord)
