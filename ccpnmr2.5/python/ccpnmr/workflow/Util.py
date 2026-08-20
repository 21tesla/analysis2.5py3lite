"""

In this file: high level stuff to handle workflow, project input.

Ask for code from partners to compare output for each component.

Need setup for CODE EXECUTION!! Similar to PDBe analysis stuff

MOSTLY also CCPN object tracking, and I/O

Real workflow examples are in python/testCode/ccpnmr/workflow

ALSO have list of components that can be tested in each software?
     Or something so that don't try to run something if it doesn't make sense?

"""

import time, tarfile, gzip, os, shutil

from subprocess import Popen, PIPE

from memops.api import Implementation

from memops.general.Io import loadProject, saveProject, packageProject

from ccp.general.Util import setCurrentStore

from ccpnmr.workflow.Constants import programList
from ccpnmr.format.general.Util import createSelection

class WorkFlow:
  
  """
  
  Prototype workflow creation wrapper for Extend-NMR.
  
  Tracks main CCPN objects, has core functions that can be used by program-specific code.
  
  """
  
  #
  # TODO Connect this to command line arguments? See eNmr.convertCasdNmrToCcpn
  #

  WorkFlowError = StandardError
  
  def __init__(self,**keywds):
  
    self.initialiseWorkflow(**keywds)

  def initialiseWorkflow(self, ccpnProject=None,          # One of these top three is obligatory!
                               ccpnProjectTgz = None,
                               identifier=None,
                               useGui=False,
                               overwrite=True,
                               initializeObjects=True,    # Set to False if don't want to automatically set NmrProject and constraint stuff. Should probably do this differently anyway! More options...
                               verbose=True):
  
    self.verbose = verbose
  
    #
    # Need either an existing CCPN project or an identifier for it - note that this
    # can be set as part of the class if required.
    #
    
    if not hasattr(self,'identifier') and identifier:
      self.identifier = identifier
    
    if not ccpnProject and not self.identifier and not ccpnProjectTgz:
      raise self.WorkFlowError("No CCPN project or identifier given - aborting workflow.")
  
    #
    # Set day/time when workflow was run
    #
    
    self.timeFlag = self.getTimeString()
    
    #
    # Set CCPN project - this is obligatory, and the CCPN project is unpacked if
    # tgz file. Is this smart? Should I always start from a loaded project?
    # And always pack it (temporary file) if have to send off? Probably best!
    #
    
    self.ccpnProject = None
    
    if ccpnProject:
      self.ccpnProject = ccpnProject
    elif ccpnProjectTgz:
      self.ccpnProject = self.unpackCcpnProjectTgz(ccpnProjectTgz)
    else:
      
      try:
        tmpCcpnProject = loadProject(self.identifier)
      except:
        tmpCcpnProject = None
      
      if tmpCcpnProject:
        if not overwrite:
          print ("  Warning: using existing project!")
          self.ccpnProject = tmpCcpnProject
        else:
          print ("  Warning: overwriting existing project!")
          shutil.rmtree(self.identifier)
      
      if not self.ccpnProject:
        self.ccpnProject = Implementation.MemopsRoot(name = self.identifier)
   
    #
    # Set other CCPN objects (if required)
    #
    
    if initializeObjects:
      # NMR project
      setCurrentStore(self.ccpnProject,'NmrProject')
      self.nmrProject = self.ccpnProject.currentNmrProject

      nmrConstraintStores = self.ccpnProject.sortedNmrConstraintStores()
      if not nmrConstraintStores:
        self.nmrConstraintStore = None
      else:
        if len(nmrConstraintStores) > 1:
          raise self.WorkFlowError("Only one NMR constraint store allowed for running workflow. Please modify input.")
        self.nmrConstraintStore = nmrConstraintStores[0]
    
    # TODO nmrEntry?
    
    #
    # Set graphical interface, if required
    #
    
    self.useGui = useGui
    if useGui:
      import Tkinter 
      self.guiRoot = Tkinter.Tk()
    else:
      self.guiRoot = None
      
    #
    # Run other <programShortCode>__init__
    #
    
    for program in programList:
      programInit = '%s__init__' % program.lower()

      if hasattr(self,programInit):
        getattr(self,programInit)()
      
  def getTimeString(self):
      
    return time.strftime("%Y-%m-%d %H:%M:%S")
        
  def getFileTimeString(self):
      
    return time.strftime("%Y.%m.%d_%H.%M")
        
  def streamData(self):
    
    #
    # So in here combine the components from Cing.py, Fc.py, ...
    #
  
    print "Please define the data stream for the workflow in a subclass component of WorkFlow!"
  
  def getResults(self):
  
    print "Get results not defined yet"


  #
  # Generally useful methods
  #
      
  def reportUnlinkedResonances(self,resonanceList):
    
    unlinkedRes = 0
    for resonance in resonanceList:
      if not resonance.resonanceSet:
        unlinkedRes += 1
        
    if unlinkedRes:
    
      if resonanceList[0].className == 'Resonance':
        parentText='nmrProject'
      else:
        parentText='nmrConstraintStore'
        
      print
      print "Warning: %d (out of %d) resonances not linked in %s." % (unlinkedRes,len(resonanceList),parentText)
      print
  
  #
  # Code to (un)pack CCPN projects into/from .tgz
  #
  
  def makeCcpnProjectTgz(self,baseNameSuffix="",outputDir=None):
  
    #
    # Set up
    #
   
    userData = self.ccpnProject.findFirstRepository(name='userData')
    ccpnProjectDir = userData.url.path
    if not outputDir:
      outputDir = ccpnProjectDir
 
    tarFileBaseName = os.path.join(outputDir + baseNameSuffix)

    tarFileName = "%s.tar" % tarFileBaseName
    tgzFileName = "%s.tgz" % tarFileBaseName
    
    #
    # Re-save before tarring, in case of changes
    #
    
    if baseNameSuffix:
      saveProject(self.ccpnProject, newPath = tarFileBaseName, removeExisting = True)
      ccpnProjectDir = userData.url.path
      
    else:
      self.ccpnProject.saveModified()
    
    # TJS: Use Wayne's packageProject, which only packages uses .xml files
       
    packageProject(self.ccpnProject, tarFileBaseName)
       
    """
    #
    # Now go to directory containing CCPN project to make .tar file, do not want full path for project!
    #
    
    curDir = os.getcwd()
    (mainCcpnDir,projectNameDir) = os.path.split(ccpnProjectDir)
    os.chdir(mainCcpnDir)

    tar = tarfile.open(tarFileName, "w")
    tar.add(projectNameDir)
    tar.close()
    
    os.chdir(curDir)
    
    #
    # Now make .tgz, remove .tar file
    #

    inFile = open(tarFileName,'rb')
    outFile = gzip.open(tgzFileName,'wb')
    outFile.writelines(inFile)
    outFile.close()
    inFile.close()
    
    os.remove(tarFileName)
    """
    
    print "  Created .tgz archive %s for CCPN project %s." % (tgzFileName,self.ccpnProject.name)

    return tgzFileName

  def unpackTgzFile(self, tgzFileName, unpackDir=None, removeFile=False, excludeFiles=None, extractFiles=None):

    """
    Unpack a .tgz or .tar.gz file with Popen and return the output from tar.
    
    Set unpackDir to the directory where the unpacked files need to go.
    
    """
        
    tgzFileName = os.path.expanduser(tgzFileName)
    
    # This has to go first.
    if unpackDir:
      extraArgs = ['-C',unpackDir]
    else:
      extraArgs = []

    if extractFiles:
      extraArgs = extraArgs + extractFiles

    if excludeFiles:
      extraArgs = extraArgs + ['--exclude'] + excludeFiles
      
    print ' '.join(['tar','xvfz',tgzFileName] + extraArgs)
    process = Popen(['tar','xvfz',tgzFileName] + extraArgs, stdin=PIPE, stdout=PIPE, close_fds=True)
    (pipe, stdin) = (process.stdout, process.stdin)

    textOutput = pipe.read()

    pipe.close()
    stdin.close()
    
    if removeFile:
      os.remove(tgzFile)

    return textOutput

  def unpackCcpnProjectTgz(self,tgzFileName,removeFile = False):
  
    """
    Unpack a CCPN project and read it
    """
  
    textOutput = self.unpackTgzFile(tgzFileName)

    lines = textOutput.split("\n")
    ccpnProjectId = lines[0]
      
    ccpnProject = loadProject(ccpnProjectId)
      
    return ccpnProject

  def unpackZipFile(self, zipFileName, unpackDir=None, removeFile=False, excludeFiles=None, extractFiles=None):

    """
    Unpack a .zip file with Popen and return the output from tar.
    
    Set unpackDir to the directory where the unpacked files need to go.
    """
        
    zipFileName = os.path.expanduser(zipFileName)

    
    if extractFiles:
      extraArgs = extractFiles
    else:
      extraArgs = []

    if unpackDir:
      extraArgs = extraArgs + ['-d',unpackDir]

    if excludeFiles:
      extraArgs = extraArgs + ['-x'] + excludeFiles
    
    process = Popen(['unzip','-o',zipFileName] + extraArgs, stdin=PIPE, stdout=PIPE, close_fds=True)
    (pipe, stdin) = (process.stdout, process.stdin)

    textOutput = pipe.read()

    pipe.close()
    stdin.close()
    
    if removeFile:
      os.remove(zipFile)

    return textOutput

  #
  # Code to tag objects - will become obsolete with NmrCalc package introduction
  #

  def tagCcpnObjectForExternalExec(self,ccpnObject,softwareName):
  
    tagCcpnObjectForExternalExec(ccpnObject,softwareName,self.getFileTimeString())
    
    if self.verbose:
      (selectionList,selectionDict) = createSelection([ccpnObject])
      print "  Tagging %s object %s" % (ccpnObject.className,selectionList[0])
      
  def tagAll(self,softwareName):
  
    #
    # Tag all chains, constraintLists, peakLists and shiftLists
    # Only relevant for automatic execution, simple projects, and testing
    #
    
    for molSystem in self.ccpnProject.sortedMolSystems():
      for chain in molSystem.sortedChains():
        self.tagCcpnObjectForExternalExec(chain,softwareName)
    
    for constraintList in self.nmrConstraintStore.sortedConstraintLists():
      self.tagCcpnObjectForExternalExec(constraintList,softwareName)
    
    for nmrExp in self.nmrProject.sortedExperiments():
      
      isDistance = False
      for expTransfer in nmrExp.expTransfers:
        if expTransfer.transferType == 'through-space':
          isDistance = True
          break
      
      if not isDistance and nmrExp.refExperiment:
        for expGraph in nmrExp.refExperiment.nmrExpPrototype.expGraphs:
          for expTransfer in expGraph.expTransfers:
            if expTransfer.transferType == 'through-space':
              isDistance = True
              break
          if isDistance:
            break
      
      if isDistance:
        for dataSource in nmrExp.sortedDataSources():
          for peakList in dataSource.sortedPeakLists():
            if peakList.peaks:
              self.tagCcpnObjectForExternalExec(peakList,softwareName)
          
    for shiftList in self.nmrProject.findAllMeasurementLists(className = 'ShiftList'):
      self.tagCcpnObjectForExternalExec(shiftList,softwareName)

#
# Generic functions - should go elsewhere eventually
#

from ccp.general.Util import setUniqueAppData

executedKw = 'executed_'

def tagCcpnObjectForExternalExec(ccpnObject,softwareName,dateTag):
  
  #
  # So the idea is that this tag will be set to True on execution of the software!
  #
  
  setUniqueAppData('AppDataBoolean',ccpnObject,softwareName,'%s%s' % (executedKw,dateTag), False)

def findCcpnObjectForExternalExec(ccpnObject,softwareName):

  date = None

  appData = ccpnObject.findAllApplicationData(application=softwareName)
  for appDatum in appData:
    if appDatum.keyword.count(executedKw):
      # Only look for tags that are not yet 'executed'
      if appDatum.value == False:
        date = appDatum.keyword.replace(executedKw,'')
        break
  
  return date
        
if __name__ == '__main__':

  workFlow = WorkFlow(identifier="test",useGui=True)
  workFlow.streamData()
  workFlow.getResults() # Is this necessary?  
  
