
from memops.universal.Io import joinPath

from ccpnmr.workflow.Constants import importDataDir
from ccpnmr.workflow.Util import WorkFlow

from ccpnmr.format.general.Conversion import FormatConversion

from ccpnmr.format.process.sequenceCompare import SequenceCompare

class FcWorkFlow(WorkFlow):

  FcWorkFlowError = StandardError

  #
  # componentList contains the data elements that this particular workflow bit can handle as input.
  # If empty, means it can handle anything.
  #
  # TODO: should I split this up into import/export?
  #
  
  componentList = []

  #
  # These are class-specific linkResonances variables that can be set for a specific
  # import.
  #
  
  specificResNameMappings = {}
  forceChainMappings = {}
  useCommonNames = False
  useIupacMatching = False
  
  def fc__init__(self):
            
    """
    Create formatConversion object and set NMR constraint store
    """
    
    self.formatConversion =  FormatConversion(ccpnProject = self.ccpnProject, guiRoot = self.guiRoot)
    
    if not hasattr(self,'nmrConstraintStore'):
      self.nmrConstraintStore = None


  def fcImportValidateAllData(self):
  
    """
    Import the data for a particular project - this is customised in subclass,
    and validate the information
    """
    
    self.fcImportAllData()
    
    #
    # Make sure nmrConstraintStore is set - TODO Where does this belong logically?
    #
    # TODO use formatConversion objects? In a way has code parallel to this...
    #

    if not self.nmrConstraintStore and self.ccpnProject.currentNmrConstraintStore:
      self.nmrConstraintStore = self.ccpnProject.currentNmrConstraintStore
    
    #
    # Validate if data OK.
    #
    
    self.fcValidateProjectData()
  
  def fcImportAllData(self):

    print "Please define an importProjectData function in a subclass component of FcWorkFlow!"
    
  def fcImportFile(self,*args,**keywds):
  
    """
    Short method name wrapper for file import
    """
  
    self.formatConversion.importFile(*args,**keywds)
    
  def fcGetFormatNameSuggestion(self,informationType,filePath):

    """
    Code to get suggestions for the file format. Will exit when done as default behaviour
    so the correct format can be manually set for this file.
    """
  
    formatNameSuggestions = self.formatConversion.determineFormatNamesForFile(informationType,filePath)

    print "Format suggestions for file %s containing %s information:" % (filePath,informationType)
    print formatNameSuggestions
    # TODO this behaviour should probably be settable
    exit(-1)
    
  def fcExportFile(self,*args,**keywds):

    """
    Short method name wrapper for file export
    """
  
    self.formatConversion.exportFile(*args,**keywds)

  def fcValidateProjectData(self):
  
    """
    Check for unlinked resonances, and report
    """
    
    self.reportUnlinkedResonances(self.nmrProject.sortedResonances())
    
    if self.nmrConstraintStore:
      self.reportUnlinkedResonances(self.nmrConstraintStore.sortedFixedResonances())


  def fcConnectResonancesToAtoms(self,importFormatName,chains):
  
    """
    Run linkResonances for an import format name and a defined set of chains
    """

    sequenceComparison = SequenceCompare()

    #
    # Set CCPN info
    #

    sequenceComparison.createCcpnChainInformation(chains)

    #
    # Get info from the resonances - NOTE this should be run after every individual import!!
    #

    resonances = self.formatConversion.getFormatClass(importFormatName).newResonances
    
    if not resonances:
      print "\n  Warning: No new resonances created during this import, skipping resonance connecting.\n"
      return None

    sequenceComparison.getFormatFileInformation(resonances,importFormatName)

    if not sequenceComparison.formatFileResidueDict.keys():
      print "\n  Warning: No format chain information available, skipping resonance connecting.\n"
      return None

    sequenceComparison.createFormatFileChainInformation()

    #
    # Now run the comparison...
    #

    if not self.forceChainMappings:
      self.forceChainMappings = sequenceComparison.compareFormatFileToCcpnInfo()

      print "\n*** Chain mappings set by alignment information ***\n"
    else:
      print "\n*** Chain mappings defined by user ***\n"

    print self.forceChainMappings

    #
    # Reset list of resonances for FormatClass!
    #

    self.formatConversion.getFormatClass(importFormatName).newResonances = []

    #
    # Now run linkResonances with automatic mapping
    #

    self.formatConversion.linkResonances(

        forceChainMappings=       self.forceChainMappings,
        useLinkResonancePopup=    True,
        setSingleProchiral=       True,                    # Assume that, e.g. HB2 is not same as HB3 if only info for HB3.
        setSinglePossEquiv=       False,                   # Assume that, e.g. HD2 is same as HD1 for Phe, Tyr if only info for HD1.
        specificResNameMappings=  self.specificResNameMappings,
        useCommonNames=           self.useCommonNames,
        useIupacMatching=         self.useIupacMatching
        )

    return self.forceChainMappings

  def fcSetPeaksInformation(self,formatName,filePath,addKeywords):
  
    """
    Handling of peak lists - can use specific addKeywords that make it easier to select
    reference experiments and set the peak list to reference experiment dimension mappings.
    """

    from ccpnmr.format.general.Util import getRefExpFromOldExpType

    #
    # Step 1, preparse the file and get info out
    #
    
    (fileRead,fileInformation) = self.formatConversion.preparseFile('peaks', formatName, filePath, addKeywords=addKeywords) 

    if not fileRead:
      raise self.FcWorkFlowError("Could not read peak list file %s in format %s:\n%s" % (filePath,formatName,fileRead))
    else:
      peakListNames = fileInformation.keys()
      
      if len(peakListNames) > 1:
        if not addKeywords.has_key("peakListName"):
          raise self.FcWorkFlowError("Multiple peak lists in file %s, and no specific 'peakListName' defined in addKeywords." % (filePath))
          
        else:
          peakListName = addKeywords['peakListName']
          
          if peakListName not in peakListNames:
            raise self.FcWorkFlowError("Multiple peak lists in file %s, and 'peakListName' %s not in list of peak list names from file." % (filePath,peakListName))
          
      else:
        peakListName = peakListNames[0]

      peakPpmRanges = fileInformation[peakListName]['peakPpmRanges']
      dimCodes =      fileInformation[peakListName]['dimensionCodes']
    
    #
    # Step 2, set the reference experiment.
    #
    
    if addKeywords.has_key('oldExpType'):
      refExperiment = getRefExpFromOldExpType(self.formatConversion.ccpnProject,addKeywords['oldExpType'])

      if not refExperiment:
        raise self.FcWorkFlowError("Could not find reference experiment 'old' name %s." % (addKeywords['oldExpType']))

      del(addKeywords['oldExpType'])

    elif addKeywords.has_key('expTypeInfo'):
      (expPrototypeName,refExpName) = addKeywords['expTypeInfo']

      refExperiment = None
      
      expPrototype = project.findFirstNmrExpPrototype(name = expPrototypeName)    
      if expPrototype:
        refExperiment = expPrototype.findFirstRefExperiment(name = refExpName)
    
      if not refExperiment:
        raise self.FcWorkFlowError("Could not find reference experiment with (prototype,experiment) names %s." % (addKeywords['expTypeInfo']))
       
      del(addKeywords['expTypeInfo'])
  
    else:
      raise self.FcWorkFlowError("No reference experiment defined for %s:\nNeed to set an 'oldExpType' or 'expTypeInfo' variable in addKeywords in self.dataFiles." % (filePath))

  
    addKeywords['refExperiment'] = refExperiment

    refExperimentInfo = self.formatConversion.getRefExperimentInfo(addKeywords['refExperiment'])

    #
    # Step 3 set the peak list to reference experiment dimension mapping, if necessary
    # 
    # uniqueMatches gives the unique (shift matches) connection between the CCPN dimension (as number) and the peak dimension (as list index)
    # matches       gives all possible shift matches
    #
    # Basically, web side has to allow all possibilities in 'matches', but should set the 'uniqueMatches' match as default.   
    #
    
    # Only run this if not manually set!

    if not addKeywords.has_key('expDimToPeakDim') or not addKeywords['expDimToPeakDim']:
    
      (matches,uniqueMatches) = self.formatConversion.matchFormatAndCcpnPeakDims(peakPpmRanges,dimCodes,refExperimentInfo['ppmRange'])

      if uniqueMatches:
        addKeywords['expDimToPeakDim'] = uniqueMatches
      else:
        raise self.FcWorkFlowError("No reference experiment dimension to peak list dimension mapping set or a unique match found. Possibilities are:\n%s\n\nReference experiment info is:\n%s\n\nPeak ppm ranges:\n%s\n" % (str(matches),refExperimentInfo,str(peakPpmRanges)))
  
#
# Below examples for importing data from format files.
#
# Alternatively, can set it up so you have to define inputs in the actual workflow... fits
# better with text-driven file strategy so people can define their own workflow.
#
# Can also keep both strategies in of course - maybe makes most sense.
#

class FcProject_1bc6(FcWorkFlow):

  """
  Defines a preset procedure to import external files into CCPN
  """

  importDir = joinPath(importDataDir,'test')

  def fcImportAllData(self):
    
    # These files have to be correct, no preparsing necessary.
    self.fcImportFile('sequence','dyana',joinPath(self.importDir,"1bc6.seq"))
    self.fcImportFile('distanceConstraints','dyana',joinPath(self.importDir,"1bc6.restr"))
    self.fcImportFile('coordinates','pdb',joinPath(self.importDir,"1bc6.pdb"), addKeywords = {'forceDefaultChainMapping': True})

    self.formatConversion.linkResonances(setSingleProchiral = True)

    # This should be modifiable...
    self.ccpnProject.saveModified()
