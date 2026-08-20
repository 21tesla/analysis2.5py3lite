""" Wim Vranken 2013. Converts CASD input data to CCPN projects.
"""


import glob, os, sys, json, shutil, traceback, tempfile

from ccpnmr.workflow.Fc import FcWorkFlow
from ccpnmr.format.general.scriptHandling import ScriptHandler
from pdbe.adatah.Io import getReferenceTextFileFromHttp, getDataFromHttp

from nijmegen.CASD import Constants as casdConstants
from nijmegen.CASD import Util as casdUtil

from memops.general import Io as genIo
allDataDir = casdConstants.allDataDir
#casdInputDir = casdConstants.casdInputDir
#casdDownloadDir = casdConstants.casdDownloadDir
#casdInputPdbDir = casdConstants.casdInputPdbDir

eNmrUrl = "http://www.wenmr.eu/wenmr/"
casdNmrDataUrl = os.path.join(eNmrUrl,"casd-nmr-data-sets")

#hrefPatt = re.compile('\<a href\=\"([^\"]+)\"')           # "<a href="anythingbutadoublequote"
#hrefNamePatt = re.compile('[^ ]\"\>([^\>]+)\<(\/a|span)') # nonspace">anythingbutendangle+</a or <span
#pdbCodePatt = re.compile('\>\s*([A-Za-z0-9]{4})\s*\<\/a') # >whitespace*4alphanums whitespace*</a     
  
#contentPattern=re.compile("[^<](?:<.*?>\s*)*([^<]*)")     # anything(<something>whitespace)* capturethis<


def getCasdNmrProjectInfo(casdNmrRefFile=None):

  """
  Code to get list of CASD-NMR projects info
  """
  
  result = []

  # This file is customisable!
  if not casdNmrRefFile:
    casdNmrRefFile = os.path.join(allDataDir,'dataPage.html')
    
  
  # Get the web page...
  text = ''.join(getReferenceTextFileFromHttp(casdNmrDataUrl,casdNmrRefFile,
                   refText = "CASD-NMR data", isGzipped = False))
  
  table, links = casdUtil.parseHtmlTable(text)
  
  tags = table[0]
  for ii in range(1,len(table)):
    dd = {}
    result.append(dd)
    ll = table[ii]
    for jj,val in enumerate(table[ii]):
      dd[tags[jj]] = val
    
    dd['DataLink'] = links[ii][0]
  
  #
  return result


class ConvertCasdNmrToCcpn(ScriptHandler,FcWorkFlow):
  

  class ConvertCasdNmrError(StandardError):
    pass
  
  #
  # Information for untarring of projects, needs ability to tweak to avoid unpacking large (spectrum) files.
  #
  
  extractArchiveFiles = []
  excludeArchiveFiles = []
  
  #
  # PDB code, if available
  #
  
  identifier = None
  pdbCode = None
  pdbChainMapping = []
  
  #
  # Fixing files on the fly sometimes, use this suffix to identify them
  #
  
  fixedSuffix = 'fixed'
  
  #
  # Options for linkResonances, might need more fine-tuning here because problems might occur
  # for specific formats or even files (think NRG distance/dihedral problems!).
  #

  specificResNameMappings = {}
  forceChainMappings = {}
  useCommonNames = False
  useIupacMatching = False
  

  #
  # Information for command line options handling
  #
  
  programDescription = "Conversion of CASD-NMR project data to CCPN projects."
  programVersion = '1.0'
  programUsage = 'Usage: %prog <options>'

  optionTuple = (
  
    ('notVerbose',  'n', False,     False, None,  'Disable verbose mode.'),
    ('interactive', 'i', False,     False, None,  'Switch to interactive mode.'),
    #('noExtract',   'x', False,     False, None,  'Do not extract original archive file.'),
    #('saveCcpn',    's', False,     False, None,  'Save CCPN project in default location.'),
    ('saveCcpn',    's', True,     False, None,  'Save CCPN project in default location.'),
    
  )  
  
  def initialise(self,keywds):
    
    
    for tag in ('identifier',):
      if tag in keywds:
        setattr(self, tag, keywds.pop(tag))
      else:
        setattr(self, tag, None)
        
        
      
    # Set the identifier for the instance
    if self.identifier is None:
      raise self.ConvertCasdNmrError("No identifier passed in to %s." % self.__class__.__name__)
      #self.identifier = self.__class__.__name__
    
    self.dataDir = os.path.join(allDataDir, self.identifier[1:3], self.identifier)
    if not os.path.isdir(self.dataDir):
      os.makedirs(self.dataDir)
    
    print "Converting CASD-NMR project %s..." % self.identifier
    
    # Create log file(s)
    self.createLogFiles()

    # These are set from script handling, but do not override input settings if available
    if self.notVerbose and not keywds.has_key('verbose'):
      keywds['verbose'] = not(self.notVerbose)

    if self.interactive and not keywds.has_key('useGui'):
      keywds['useGui'] = self.interactive
        
    # Now initialise Workflow script...
    self.initialiseWorkflow(**keywds)

  def fcImportAllData(self):
    
    #
    # Get CASD-NMR data
    #
    print '### import', self.identifier
    
    self.tmpdir = tempfile.mkdtemp(dir=casdConstants.topTmpDir)
    try:
      self.unpackDir = os.path.join(self.tmpdir,'unpack')
      os.mkdir(self.unpackDir)
      
      self.unpackCasdTgz()
 
      self.fileModifier()
 
      self.importCasdData()
      
      # NB adds all data in relevant NmrProject to run 
      #casdUtil.makeNmrCalcRun(self.ccpnProject, task='CASD-NMR')
      casdUtil.makeNmrCalcRun(self.ccpnProject, task='CING')
 
      self.saveProject()

    finally:
      shutil.rmtree(self.tmpdir, ignore_errors=True)
    
    
    # Close log file(s)
    self.closeLogFiles()
  
  #def importPatches(self):
  #  print '### patches', self.identifier, self.patchesDir, os.path.exists(self.patchesDir)
  #  if os.path.exists(self.patchesDir):
  #    for ff in os.listdir(self.patchesDir):
  #      path = os.path.join(self.patchesDir, ff)
  #      print '### patches', self.importDir, path
  #      shutil.copy(path, self.importDir)
  
  
  def importPdbFileInfo(self):
  
    if self.pdbCode:
    
    
      pdbFileName = self.getInputFile('structures', ignoreErrors=True)
      if pdbFileName:
        addKeywords = {}
        if self.pdbChainMapping:
          addKeywords['forceChainMappings'] = self.pdbChainMapping

        self.fcImportFile('coordinates','pdb', pdbFileName, addKeywords=addKeywords)
      
      #from pdbe.adatah.Pdb import getPdbCode, pdbDataDir
      
      # Note: can fail if not yet released!   
      #if getPdbCode(self.pdbCode):

      #  pdbFileName = os.path.join(pdbDataDir,"%s.pdb" % self.pdbCode.lower())

      #  addKeywords = {}
      #  if self.pdbChainMapping:
      #    addKeywords['forceChainMappings'] = self.pdbChainMapping

      #  self.fcImportFile('coordinates','pdb',pdbFileName, addKeywords=addKeywords)
            
  def fileModifier(self):
  
    # Use this to modify files that are not parseable because their format is not defined
    # or misses crucial information. No point in writing parsers for such things.
    pass
  
  def getInputFile(self, subdir, ignoreErrors=False):
    
    return casdUtil.getInputFile(self.identifier, subdir,
                                 ignoreErrors=ignoreErrors)
    
    
  def unpackCasdTgz(self):
    
    casdFile = self.getInputFile('restraints')
        
    print "  Unpacking CASD-NMR project in %s..." % self.unpackDir

    if casdFile.count(".zip"):
      textOutput = self.unpackZipFile(casdFile, unpackDir=self.unpackDir, excludeFiles=self.excludeArchiveFiles, extractFiles=self.extractArchiveFiles)
    else:
      textOutput = self.unpackTgzFile(casdFile, unpackDir=self.unpackDir, excludeFiles=self.excludeArchiveFiles, extractFiles=self.extractArchiveFiles)
    
    #if not textOutput:
    #  raise self.ConvertCasdNmrError("Tar unpacking did not work for file %s." % (casdFile))

    if textOutput.count("tar:"):
      raise self.ConvertCasdNmrError("Tar unpacking error for file %s:\n%s" % (casdFile,textOutput))

    else:
      # This is subclass modifiable in case files have to be moved/renamed, ...
      self.setImportDir(textOutput)
  
  def setImportDir(self,textOutput):
    
    self.importDir = casdUtil.getLowestSubDir(self.unpackDir)
    
  def importCasdData(self):
  
    doPdbImport = False
    
    for (informationType,formatName,fileName,addKeywords) in self.dataFiles:
    
      filePath = os.path.join(self.importDir,fileName)
      
      if not formatName:
        print "Determining format for %s file %s..." % (informationType,filePath)
        self.fcGetFormatNameSuggestion(informationType,filePath)
      
      #
      # PDB file is slightly tricky - do this AFTER project, but BEFORE sequence import (and ignore the sequence import)
      #
      # This is necessary to avoid sequence duplication, but still necessary for consistent import of NMR-STAR data.
      #
      
      if self.pdbCode:
        if informationType in ('project','sequence'):
          doPdbImport = True
          
      #
      # Get information from PDB code, if available
      #
      
      if doPdbImport and informationType != 'project':

        self.importPdbFileInfo()
        
        doPdbImport = False
        
        if informationType == 'sequence':
          continue
      
      #
      # Special handling for sequences - use self.identifier as molecule name if possible.
      #
      
      if informationType == 'sequence':
        if not addKeywords.has_key('molName'):
          addKeywords['molName'] = self.identifier
     
      #
      # Special handling for peaks
      #
      
      if informationType == 'peaks':
        self.fcSetPeaksInformation(formatName,filePath,addKeywords)
        
        if 'specName' not in addKeywords:
          # Get a reasonable name for the peak list
          peakListNameItems = fileName.split('.')
          if len(peakListNameItems) > 1:
            peakListName = '.'.join(peakListNameItems[:-1])
          else:
            peakListName = fileName
 
          addKeywords['specName'] = peakListName
       
      # TODO add options here!
      self.fcImportFile(informationType,formatName,filePath,addKeywords=addKeywords)
      
      # Run linkResonances for all non-sequence and coordinate info. Make sure that you read
      # in the sequence first!!
      if informationType not in ('sequence','coordinates'):
      
        molSystem = self.ccpnProject.findFirstMolSystem()        
        if not molSystem:
          raise self.ConvertCasdNmrError("No molecular system available - make sure to import sequence or coordinates first.")

        chains = molSystem.chains
        if not chains:
          raise self.ConvertCasdNmrError("No chains available - make sure to import sequence or coordinates first.")

        self.fcConnectResonancesToAtoms(formatName,chains=chains)

  def createLogFiles(self):
    
    logDir = os.path.join(self.dataDir, 'log_extractEntry')
    if not os.path.isdir(logDir):
      os.makedirs(logDir)
    logFilePath = os.path.join(logDir, '%s.%s.log' % (self.identifier,self.getFileTimeString()))        
    self.logFileHandle = open(logFilePath,'w')
    
    sys.stdout = self.logFileHandle
    
    #self.errorLogFilePath = os.path.join(getLogfileDir(),'%s.%s.error' % (self.logFileRoot,self.currentTimeStamp))        
    #self.errorLogFileHandle= open(self.errorLogFilePath,'w')

  def closeLogFiles(self):
  
    self.logFileHandle.close()    
    
    sys.stdout = sys.__stdout__
    
    #self.errorLogFileHandle.close()
    
    #if not self.hasErrors:
    #  os.remove(self.errorLogFilePath)
    
  def saveProject(self):
  
    if self.saveCcpn:
      
      # Must be saved (and thus also renamed) before it can be packaged
      genIo.saveProject(self.ccpnProject, 
                        newPath=os.path.join(self.dataDir, self.identifier),
                        newProjectName=self.identifier, 
                        checkValid=True, removeExisting=True)
                        
      genIo.packageProject(self.ccpnProject, 
                           os.path.join(self.dataDir, self.identifier))
 
class AR3436A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'AR3436A_blind'
  casdUnpackDir =           'AR3436A_blind' # TODO in principle could force this to be self.identifier...
  
  pdbCode = '2KJ6'
  
  dataFiles = (
  
    #("sequence",    'sparky',   "AR3436A_seq"), Not necessary, is in bmrb file.
    ("project",     'nmrStar',  "AR3436A_0309.bmrb",                     {'version': '2.1.1'}),
    ("peaks",       'sparky',   "AR3436ANC_Nnoesy_final.list",           {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'sparky',   "AR3436A_Caliphnoesy_folded_final.list", {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'sparky',   "AR3436ANC_aromnoesy_final.list",        {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )


class CGR26A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'CGR26A_blind'
  #casdUnpackDir =       not necessary, is renamed!
  
  pdbCode = '2KPT'

  extractArchiveFiles = ['casd']
  excludeArchiveFiles = ['*.tgz']
    
  dataFiles = (
  
    ("sequence",    'sparky',   "cgr26a.aa",         {}),
    ("shifts",      'nmrStar',  "cgr26a.bmrb",       {'version': '2.1.1'}),
    ("peaks",       'xeasy',    "n.peaks",           {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'xeasy',    "ali.peaks",         {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'xeasy',    "aro.peaks",         {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )


class CtR69A(ConvertCasdNmrToCcpn):

  #casdUnpackDir = None  - no need to set, is identifier!
  
  pdbCode = '2KRU'
  
  dataFiles = (
   
    ("sequence",     'xeasy',    "ctr69a.seq",              {}),                                
    ("shifts",       'nmrStar',  "ctr69a_nmrstar31.dat",    {'version': '2.1.1'}),              
    ("peaks",        'xeasy',    "n.peaks",                 {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",        'xeasy',    "cali.peaks",              {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",        'xeasy',    "caro.peaks",              {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )
    

class ET109Ared(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = "ET109A_NMR_CASD"
  casdUnpackDir = 'reduced_ET109A'
  
  pdbCode = '2KKX'
  
  extractArchiveFiles = ['reduced_ET109A']
  
  dataFiles = (
  
    ("sequence",    'sparky',   "reduced_ET109A_seq",                    {}),
    ("shifts",      'nmrStar',  "reduced_ET109A.bmrb.fixed",             {'version': '2.1.1'}),
    # Not working ("rdc",         'cyana',    "reduced_ET109A_RDC.in",                 {}),
    ("peaks",       'xeasy',    "reduced_ET109A_Nnoesy.peaks",           {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'xeasy',    "reduced_ET109A_Cnoesy.peaks",           {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'xeasy',    "reduced_ET109A_Cnoesy_aro.peaks",       {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )
  
  def fileModifier(self):
   
    # The chemical shift file does not have a loop_ with info on the column content.

    origFilePath = os.path.join(self.importDir,'reduced_ET109A.bmrb')
    fixedFilePath = "%s.%s" % (origFilePath,self.fixedSuffix)

    if not os.path.exists(fixedFilePath):

      fin = open(origFilePath)
      lines = fin.readlines()
      fin.close()

      fout = open(fixedFilePath,'w')
      fout.write("""   
    loop_
      _Atom_shift_assign_ID
      _Residue_seq_code
      _Residue_label
      _Atom_name
      _Atom_type
      _Chem_shift_value
      _Chem_shift_value_error
      _Chem_shift_ambiguity_code
      
      """)

      for line in lines:
       fout.write(line)

      fout.write("\nstop_\n")


class ET109Aox(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = "ET109A_NMR_CASD"
  casdUnpackDir = 'oxidized_ET109A'

  pdbCode = '2KKY'
  
  extractArchiveFiles = ['oxidized_ET109A']
  
  dataFiles = (
  
    ("sequence",    'sparky',   "oxidized_ET109A_seq",                    {}),
    ("shifts",      'nmrStar',  "oxidized_ET109A.bmrb.fixed",             {'version': '2.1.1'}),
    # Not working ("rdc",         'cyana',    "oxidized_ET109A_RDC.in",                 {}),
    ("peaks",       'xeasy',    "oxidized_ET109A_Nnoesy.peaks",           {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'xeasy',    "oxidized_ET109A_Cnoesy.peaks",           {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'xeasy',    "oxidized_ET109A_Cnoesy_aro.peaks",       {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )
  
  def fileModifier(self):
   
    # The chemical shift file does not have a loop_ with info on the column content.

    origFilePath = os.path.join(self.importDir,'oxidized_ET109A.bmrb')
    fixedFilePath = "%s.%s" % (origFilePath,self.fixedSuffix)

    if not os.path.exists(fixedFilePath):

      fin = open(origFilePath)
      lines = fin.readlines()
      fin.close()

      fout = open(fixedFilePath,'w')
      fout.write("""   
    loop_
      _Atom_shift_assign_ID
      _Residue_seq_code
      _Residue_label
      _Atom_name
      _Atom_type
      _Chem_shift_value
      _Chem_shift_value_error
      _Chem_shift_ambiguity_code
      
      """)

      for line in lines:
       fout.write(line)

      fout.write("\nstop_\n")


class HR5537A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'HR5537A_bind'
  casdUnpackDir =           'HR5537A_bind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2KK1'
  
  dataFiles = (
  
    #("sequence",    'sparky',   "HR5537A_bind/HR5537A.seqs", {}),
    #("shifts",    'nmrStar',   "HR5537A_B800_final.bmrb", {}),
    ("project",     'nmrStar',  "16349_km_rohcoomk.str",            {'version': '2.1.1'}),
    
    ("peaks",       'xeasy',   "HR5537A_B800_simNOESY_15N_unfolded_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'xeasy',   "HR5537A_B800_simNOESY_15N_unfolded.peaks",      {'oldExpType': 'noesy_hsqc_HNH.hhn'}),

    ("peaks",       'xeasy',   "HR5537A_B800_simNOESY_13C_unfolded_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'xeasy',   "HR5537A_B800_simNOESY_13C_unfolded.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc'}),

    ("peaks",       'xeasy',   "HR5537A_B800_D2ONOESY_13Caliph_unfolded_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'xeasy',   "HR5537A_B800_D2ONOESY_13Caliph_unfolded.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'xeasy',   "HR5537A_B800_simNOESY_13Caliph_unfolded.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )


class NeR103A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'NeR103A_blind'
  casdUnpackDir =           'NeR103A_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2KPM'

  excludeArchiveFiles = ['*.tar']

  dataFiles = (
  
    ("project",     'nmrStar',  "16560_NeR103A_3.1_updated.bmrb",      {'version': '2.1.1'}),
    ("peaks",       'sparky',   "NeR103A_CASD_n_noesy.list",          {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'sparky',   "NeR103A_CASD_ali_noesy.list",        {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'sparky',   "NeR103A_CASD_aro_noesy.list",        {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )


class PGR122A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'PGR122A_blind'
  casdUnpackDir =           'PGR122A_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2KMM'

  dataFiles = (
  
    ("sequence",    'fasta',    "sequence.fasta.fixed",        {}),
    ("shifts",      'nmrStar',  "pgr122_bmrb_3.1.txt.fixed",   {'version': '3.1'}),
    ("peaks",       'sparky',   "Nnoesy.list",                 {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'sparky',   "Cnoesy_aliph_850MHz.list",    {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'sparky',   "Cnoesy_arom.list",            {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )


  def fileModifier(self):
  
    # No sequence file, creating this

    fixedFilePath = os.path.join(self.importDir,'sequence.fasta.%s' % self.fixedSuffix)
  
    if not os.path.exists(fixedFilePath):
      fout = open(fixedFilePath,'w')
      fout.write("> File written by convertCasdNmrToCcpn because not present.\n")
      fout.write("MEVMVFTPKGEIKRLPQGATALDFAYSLHSDLGDHCIGAKVNHKLVPLSYVLNSGDQVEVLSSKSLEHHHHHH\n")
      fout.close()
   
    # The chemical shift file contains nonsense.

    origFilePath = os.path.join(self.importDir,'pgr122_bmrb_3.1.txt')
    fixedFilePath = "%s.%s" % (origFilePath,self.fixedSuffix)

    if not os.path.exists(fixedFilePath):

      fin = open(origFilePath)
      lines = fin.readlines()
      fin.close()
      
      # Crap lines in file, remove them.
      newLines = []
      for line in lines:
        cols = line.split()
        if line[0] != '#' and len(cols) > 6:
          serial = int(cols[0])
          if serial >= 810:
            continue
        
        newLines.append(line)
      
      # Write out fixed version
      fout = open(fixedFilePath,'w')
      for line in newLines:
        fout.write(line)
      fout.close()


class VpR247(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'VpR247_blind'
  casdUnpackDir =           'VpR247_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2KIF'
  
  dataFiles = (
  
    ("sequence",    'xeasy',    "VpR247.seq",                                 {}),
    ("shifts",      'nmrStar',  "VpR247_B800_final.bmrb.fixed",               {'version': '3.1'}),
    ("peaks",       'sparky',   "VpR247_B800_simNOESY_15N_refined.list",      {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'sparky',   "VpR247_B800_simNOESY_13Caliph_refined.list", {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'sparky',   "VpR247_B800_aromNOESY_refined.list",         {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )


  def fileModifier(self):
   
    # This one mixes 3.1 info with a 2.1 shift info list, and has inappropriate comments in the loop_

    origFilePath = os.path.join(self.importDir,'VpR247_B800_final.bmrb')
    fixedFilePath = "%s.%s" % (origFilePath,self.fixedSuffix)

    if not os.path.exists(fixedFilePath):

      fin = open(origFilePath)
      lines = fin.readlines()
      fin.close()

      fout = open(fixedFilePath,'w')
      fout.write("""   
    loop_
      _Atom_shift_assign_ID
      _Residue_seq_code
      _Residue_label
      _Atom_name
      _Atom_type
      _Chem_shift_value
      _Chem_shift_value_error
      _Chem_shift_ambiguity_code
      
      """)

      shiftInfoStarted = False
      for line in lines:
        if line.count("# 03-18-09; 2.1 format"):
          shiftInfoStarted = True
        
        if shiftInfoStarted:
          if line.count("#"):
            line = line[:line.index("#")]
        
          fout.write(line)

      fout.write("\nstop_\n")


class AtT13(ConvertCasdNmrToCcpn):

  pdbCode = '2KNR'
  pdbChainMapping = [('A', 'A', 1, -3)]
  
  dataFiles = (
   
   # ("sequence",     'xeasy',    "proteinSequence.dat",     {}),                                
    ("project",      'nmrStar',  "atc0905.nmrstr",          {'version': '2.1.1'}),              
    ("peaks",        'xeasy',    "N15edited.peaks",         {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",        'xeasy',    "C13edited_h2o.peaks",     {'oldExpType': 'noesy_hsqc_HCH.hhc'}),        
    ("peaks",        'xeasy',    "C13_aromatic.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc'}),        
      
  )


class HR6470A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'HR6470A_blind'
  casdUnpackDir =           'HR6470A_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2L9R'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    ("peaks",       'xeasy',   "nnoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "cnoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "anoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )


class HR6430A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'HR6430A_blind'
  casdUnpackDir =           'HR6430A_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2LA6'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    ("peaks",       'xeasy',   "nnoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "cnoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "anoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )

    
class HR5460A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'HR5460A_blind'
  casdUnpackDir =           'HR5460A_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2LAH'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    ("peaks",       'xeasy',   "nnoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "cnoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "anoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )


class OR36(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'OR36_blind'
  casdUnpackDir =           'OR36_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2LCI'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    ("peaks",       'xeasy',   "nnoe_raw.peaks",    {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "alinoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "aronoe_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )


class OR135(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'OR135_blind'
  casdUnpackDir =           'OR135_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2LN3'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    # NOTE: to get mappings, just import peak lists manually with FC to be sure it's correct...
    ("peaks",       'xeasy',   "simnoe_15N_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "simnoe_13C_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "aronoe_raw.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )



class StT322(ConvertCasdNmrToCcpn):

  #casdUnpackDir = None  - no need to set, is identifier!
  
  pdbCode = '2LOJ'

  casdArchiveFileNameBase = 'StT322_blind'
  casdUnpackDir =           'StT322_blind' # TODO in principle could force this to be self.identifier...
  
  dataFiles = (
  
    ("sequence",    'fasta',    "sequence.fasta.fixed",           {}),   
    ("shifts",      'nmrStar',  "StT322_BMRB_18214.str.fixed",    {'version': '3.1'}),# Original hand-generated, not real full NMR-STAR file!
    ("peaks",       'sparky',   "autopeak_15Nnoe.list",           {'oldExpType': 'noesy_hsqc_HNH.hhn'}),
    ("peaks",       'sparky',   "autopeak_13Cnoe_h2o.list",       {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'sparky',   "autopeak_13Cnoe_d2o.list",       {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
    ("peaks",       'sparky',   "autopeak_13Caronoe.list",        {'oldExpType': 'noesy_hsqc_HCH.hhc'}),
      
  )



  def fileModifier(self):
  
    # No sequence file, creating this

    fixedFilePath = os.path.join(self.importDir,'sequence.fasta.%s' % self.fixedSuffix)
  
    if not os.path.exists(fixedFilePath):
      fout = open(fixedFilePath,'w')
      fout.write("> File written by convertCasdNmrToCcpn because not present.\n")
      fout.write("MSRMDNTELPHPKEIDNETLLPAAERRVNSQALLGPDGKVIIDHNGQEYLLRKTQAGKLLLTK\n")
      fout.close()
   
    # The chemical shift file contains nonsense.

    origFilePath = os.path.join(self.importDir,'StT322_BMRB_18214.str')
    fixedFilePath = "%s.%s" % (origFilePath,self.fixedSuffix)

    if not os.path.exists(fixedFilePath):

      fin = open(origFilePath)
      lines = fin.readlines()
      fin.close()
      
      # info missing, add it
      newLines = []
      newLines.append("data_title\n")
      newLines.append("save_assigned_chemical_shifts\n")
      newLines.append("   _Assigned_chem_shift_list.Sf_category              assigned_chemical_shifts\n")
      newLines.append("   _Assigned_chem_shift_list.Entry_ID      1\n")
      newLines.append("   _Assigned_chem_shift_list.ID            1\n")

      for line in lines:
        if line.count("Comp_ID"):
          newLines.append("        _Atom_chem_shift.Comp_index_ID\n")          
        newLines.append(line)
      newLines.append("save_\n")
      
      # Write out fixed version
      fout = open(fixedFilePath,'w')
      for line in newLines:
        fout.write(line)
      fout.close()

class HR2876C(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'HR2876C_2m5o_casd_nmr'
  casdUnpackDir =           'HR2876C_2m5o' # TODO in principle could force this to be self.identifier...

  pdbCode = '2M5O'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    # NOTE: to get mappings, just import peak lists manually with FC to be sure it's correct...
    # NB the following are taken from HR2867B, assuming they are the same.
    ("peaks",       'xeasy',   "simnoe_15N_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "simnoe_13C_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "aronoe_raw.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )


class HR2876B(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'HR2876B_blind'
  casdUnpackDir =           'HR2876B_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2LTM'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    # NOTE: to get mappings, just import peak lists manually with FC to be sure it's correct...
    ("peaks",       'xeasy',   "simnoe_15N_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "simnoe_13C_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "aronoe_raw.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    #("peaks",       'xeasy',   "simnoe_15N_final.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    #("peaks",       'xeasy',   "simnoe_13C_final.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    #("peaks",       'xeasy',   "aronoe_final.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    #("rdcConstraints", 'cyana',   "cya.rdc",      {}),
      
  )


class YR313A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'YR313A_blind'
  casdUnpackDir =           'YR313A_blind' # TODO in principle could force this to be self.identifier...

  pdbCode = '2LTL'
  
  dataFiles = (
  
    ("project",     'nmrStar',  "bmrb21.str", {'version': '2.1.1'}),
    
    # NOTE: to get mappings, just import peak lists manually with FC to be sure it's correct...
    ("peaks",       'xeasy',   "simnoe_15N_raw.peaks",  {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "simnoe_13C_raw.peaks",  {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "aronoe_raw.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )

    
class MiR12(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'MiR12_for_CASD-NMR'
  casdUnpackDir =           'MiR12_for_CASD-NMR' # TODO in principle could force this to be self.identifier...

  pdbCode = None
  
  dataFiles = (
  
    ("project",     'nmrStar',  "MiR12_bmrb21.str", {'version': '2.1.1'}), # Modified file!
    
    ("peaks",       'xeasy',   "RAW_peak_lists/RAW_n15_xe.peaks",      {'oldExpType': 'noesy_hsqc_HNH.hhn'}),# 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "RAW_peak_lists/RAW_ali_xe.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "RAW_peak_lists/RAW_aro_xe.peaks",      {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}})
      
  )

class MiR12_final(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'MiR12_for_CASD-NMR'
  casdUnpackDir =           'MiR12_for_CASD-NMR' # TODO in principle could force this to be self.identifier...

  pdbCode = None
  
  dataFiles = (
  
    ("project",     'nmrStar',  "MiR12_bmrb21.str", {'version': '2.1.1'}), # Modified file!
    
    ("peaks",       'xeasy',   "final_peak_lists/n15_xe.peaks",   {'oldExpType': 'noesy_hsqc_HNH.hhn'}),# 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "final_peak_lists/ali_xe.peaks",   {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "final_peak_lists/aro_xe.peaks",   {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 2, 2: 1, 3: 3}}),
    ("peaks",       'xeasy',   "final_peak_lists/4d_xe.peaks",    {'oldExpType': 'noesy_hsqc_HCCH.hhcc', 'expDimToPeakDim': {1: 4, 2: 1, 3: 2, 4: 3}}),
      
  )


class HR8254A(ConvertCasdNmrToCcpn):

  casdArchiveFileNameBase = 'hr8254a_2m2e_casd_nmr'
  casdUnpackDir =           'hr8254a_2m2e' # TODO in principle could force this to be self.identifier...

  pdbCode = '2M2E'
  
  dataFiles = (

    ("sequence",    'fasta',    "hr8254a_fixed.seq",           {}),   
    ("shifts",      'nmrStar',  "hr8254a_bmrb_18909_mod.str",    {'version': '3.1'}),# Original hand-generated, not real full NMR-STAR file!
    
    ("peaks",       'sparky',   "hr8254a_N15noe_raw.list",    {'oldExpType': 'noesy_hsqc_HNH.hhn', 'expDimToPeakDim': {1: 1, 2: 2, 3: 3}}),
    ("peaks",       'sparky',   "hr8254a_C13noe_raw.list",    {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 1, 2: 2, 3: 3}}),
    ("peaks",       'sparky',   "hr8254a_AromNoe_raw.list",   {'oldExpType': 'noesy_hsqc_HCH.hhc', 'expDimToPeakDim': {1: 1, 2: 2, 3: 3}}),
      
  )

if __name__ == '__main__':

  import sys
  
  # get json file name
  calcDataFile = os.path.join(allDataDir, 'calcData.json') 
  calcData = json.load(open(calcDataFile))
  
  #from pdbe.adatah.CasdNmr import getCasdNmrProjectInfo
  #inputData = getCasdNmrProjectInfo()
  #json.dump(projectInfo, open(jsonFile, 'w'), sort_keys=True, indent=4)
  
  entries = ('2m5o_Org',)
  #entries = ()
  
  if not entries:
    entries = [tag for tag,val in calcData.items() if val.get('isOriginal')]
  
  for entryName in entries:
    
    targetName = calcData[entryName]['Target']
      
    #
    # Do the conversion...
    #
    converter = vars()[targetName](identifier=entryName)
    try:
      converter.fcImportAllData()
    except:
      print 'ERROR importing %s' % entryName
      traceback.print_exc(file=sys.stdout)
  
  
  #for info in inputData:
   
  #  if not (info.get('Invalid') or (targets and info['Target'] not in targets)):
 
  #    entryName = casdUtil.getEntryName(info, isOriginal=True)
  #    targetName = info['Target']
 
      
  #    #
  #    # Do the conversion...
  #    #
 
  #    converter = vars()[targetName](identifier=entryName)
 
  #    try:
  #      converter.fcImportAllData()
  #    except:
  #      print 'ERROR importing %s' % targetName
  #      traceback.print_exc(file=sys.stdout)
