""" 
"""

import os, sys, json, re, time
import tarfile, zipfile, bz2, gzip, shutil
import traceback

from ccp.general import Io as genIo
#from memops.general import Io as genIo
from nijmegen.CASD import Constants as casdConstants

allDataDir = casdConstants.allDataDir

class Extractor():
  def __init__(self, file):
    self.file = file
  
  def extractall(self, dest):
    if os.path.isdir(dest):
      destFile = os.path.splitext(os.path.basename(self.file))[0]
      destFile = os.path.join(dest, destFile)
    else:
      destFile = dest
    open(destFile, 'w').write(self.read())

class BZ2Extractor(Extractor):
  def read(self):
    return bz2.decompress(open(self.file).read())

class GZExtractor(Extractor):
  def read(self):
    return gzip.open(self.file, 'rb').read()
  

extractorClasses = (
 ('.tar.bz2', tarfile.open),
 ('.tgz', tarfile.open),
 ('.tar.gz', tarfile.open),
 ('.tar', tarfile.open),
 ('.zip', zipfile.ZipFile),
 ('.gz', GZExtractor),
 ('.bz2', BZ2Extractor),
)

  
def getEntryName(info, isOriginal=False):
  """ Get entry name from info dictionary. 
  if isOriginal get name of CASD input (...._Org)
  NBNB copy of function in cing/python/cing/NRG/CasdScripts.py
  NBNB must be consolidated
  """
  mainId = info.get('PDBcode')
  if mainId:
    mainId = mainId.lower()
  else:
    mainId = info['Target']
  
  if isOriginal or 'Submitted on' not in info:
    #Original data - put in Org
    return mainId + '_Org'
  
  else:
    entryID = info.get('EntryID', '???')
    group = info.get('Group', '????')
    return '%s_%s_%s' % (mainId, group, entryID, )
  


def parseHtmlKeyValTables(text):
  """ Parse text with a series of two-column keyword-value tables
  And put result in JSON-like structure.
  
  # NB this could be done neater with regular expressions
  """
  
  skipKeys = ('Edit Validation URLS',)
  
  result = []
  
  chunks = text.split('</table>')
  for chunk in chunks[:-1]:
    dd = {}
    rows = chunk.split('</tr>')
    for row in rows[:-1]:
      fields = row.split('<td')[1:]
      for ii in 0,1:
        field = fields[ii]
        field = field.split('</')[0]
        fields[ii] = field[field.rfind('>')+1:]
      
      # Convert text to proper value
      try:
        val = int(fields[1])
      except:
        val = None
      
      if val is None:
        try:
          val = float(fields[1])
        except:
          val = fields[1] or None
        
      key = fields[0]
      if key not in skipKeys:
        dd[key] = val
      
    result.append(dd)
  #
  return result

def parseHtmlTable(text):
  """ Find first table in text and return 1) list of contents, row by row,
  2) list of links, row by row.
  """
  #hrefPatt = re.compile('\<a href\=\"([^\"]+)\"')           # "<a href="anythingbutadoublequote"
  linkPattern = re.compile('\<a href\=\"(.+?)\"')           # "<a href="anythingbutadoublequote"
  #tablePattern = re.compile('<table.*?>(.+?)</table>', re.M)
  tablePattern = re.compile('<table.*?>(.+?)</table>')
  rowPattern = re.compile('<tr.*?>(.+?)</tr>')
  fieldPattern = re.compile('<td.*?>(.+?)</td>')
  
  # COntent pattern. NBNB does NOT strip final whitespace
  contentPattern=re.compile("\s*(?:<.*?>\s*)*([^<]*)")     # anything(<something>whitespace)* capturethis<
  
  text = ''.join(text.splitlines())
  #print text
  
  # 
  result = []
  links = []
  table = tablePattern.search(text).group()
  rows = rowPattern.findall(table)
  for row in rows:
    ll = []
    ll2 = []
    result.append(ll)
    links.append(ll2)
    row = row.replace('&nbsp;',' ')
    fields = fieldPattern.findall(row)
    for field in fields:
      x = contentPattern.match(field)
      if x:
        content = x.group(1).strip() or None
      else:
        content = None
      ll.append(content)
      x = linkPattern.search(field)
      ll2.append(x and x.group(1))
  #
  return result, links

# Too dangerous - auto clears targe directory. Do not use
#def makeClearDirectory(directory):
#  if os.path.exists(directory):
#    for ff in os.listdir(directory):
#      ff = os.path.join(directory,ff)
#      if os.path.isdir(ff):
#        shutil.rmtree(ff)
#      else:
#        os.remove(ff)
#  else:
#    os.makedirs(directory)

def extractCompressedFile(source, targetDir, dataId, okExts=()):
  """ Extract contents of compressed source into targetDir
  """
  
  # make or clear target directory
  if not os.path.exists(targetDir):
    os.makedirs(targetDir)
  
  for tag,extractor in extractorClasses:
    if source.endswith(tag):
      try:
        extractor(source).extractall(targetDir)
      except:
        print 'ERROR for %s, reading %s' % (dataId, source)
        traceback.print_exc(file=sys.stdout)
        return
      break
  else:
    for ext in okExts:
      if source.endswith(ext):
        break
    else:
      print ('### WARNING, %s Unexpected file type for: %s' %
             (dataId, os.path.basename(source)))
    shutil.copy(source, targetDir)
  

def extractToJson(htmlFile):
  """ special hack: Read an HTML file with keyword-value tables,
  and extract contents fo json file in same location
  """
  
  # get json file name
  name, ext = os.path.splitext(htmlFile)
  jsonFile = name + '.json'
  
  # get parsed values
  text = open(htmlFile).read()
  result = parseHtmlKeyValTables(text)
  
  # 
  
  refKeys = set(result[0].keys())
  for dd in result[1:]:
    if set(dd.keys()) != refKeys:
      print '~~~ Keys differ ', sorted(dd.keys())
  
  # write result to json
  json.dump(result, open(jsonFile, 'w'), sort_keys=True, indent=4)

def createDirs(dictlist):
  """ Make directories for entries, in CWD.
  """
  
  for dd in dictlist:
    target = dd['Target']
    entry = str(dd['EntryID'])
    path = os.path.join(target, entry)
    if not os.path.isdir(path):
      print '... creating ', path
      os.makedirs(path)
    print ss
  
  
def getInputFile(identifier, subdir, ignoreErrors=False):
  
  result = None 
  path = os.path.join(allDataDir, identifier[1:3], identifier,
                      '%s.input' % identifier, subdir)
  ll = os.listdir(path)
  if len(ll) == 1:
    result = os.path.join(path, ll[0])
  elif not ignoreErrors:
    raise Exception("Contents of %s not a single file: %s." % (path, ll))
  #
  return result
    

def getLowestSubDir(inputDir, followDirs=None):
  """ replace directory by subdirectory while directory contains nothing else
  If followDIrs is set, the function will go into any directory with a name on the
  list, even if the parent directory is not empty
  """
  
  result = inputDir
  
  while True:
    ll = os.listdir(result)
    if len(ll) == 1:
      name = ll[0]
    elif followDirs:
      for name in followDirs:
        if name in ll:
          break
      else:
        break
    else:
      break
        
    nextDir = os.path.join(result,name)
    if os.path.isdir(nextDir):
      result = nextDir
    else:
      break
  #
  return result
  

def replaceInSecondLines(path, src=' 1HZ  LYS', dst=' 2HZ  LYS'):
  """ Custom hack to change the second of two successive lines
  """
  
  lines = open(path).read().splitlines()
  
  for ii, line in enumerate(lines):
    if src in line:
      lines[ii+1] = lines[ii+1].replace(src, dst)
  #
  open(path,'w').write('\n'.join(lines))


def fixFourLetterHNames(path):
  """ Fixes error4 when names like 1HG1, 2HG1, 2HH1, etc. have been wrongly 
  entered as 1HG, 1HG, 1HH etc,
  """
  
  lines = open(path).read().splitlines()
  
  for ii, line in list(enumerate(lines))[:-2]:
    
    if len(line) > 50 and not line.startswith('REMARK') and line[12] in '12':
      # ATOM or HETATM records (delcared or not) with 1 or2 iposition 12
      
      # Find stretches of lines with duplicate atom identifier
      ss = line[11:28]
      nn = 1
      while lines[ii+nn][11:28] == ss:
        nn += 1
      
      # Fix the names
      if nn > 1:
        for jj in range(nn):
          ss = lines[ii+jj]
          lines[ii+jj] = ss[:12] + str(jj+1) + ss[13:15] + ss[12] + ss[16:]
  #
  open(path,'w').write('\n'.join(lines))

def moveHIndicesToSuffix(path):
  
  lines = open(path).read().splitlines()
  
  for ii, line in list(enumerate(lines))[:-2]:
    if (len(line) > 50 and not line.startswith('REMARK') 
         and line[12] in '123' and line[13] == 'H'):
    
      atName = line[12:16]
      if atName[2] == ' ':
        atName = ''.join((' ',atName[1], atName[0], ' '))
      elif atName[3] == ' ':
        atName = ''.join((' ',atName[1], atName[2], atName[0]))
      else:
        atName = atName[1:] + atName[0]
    
      lines[ii] = line[:12] + atName + line[16:]
  #
  open(path,'w').write('\n'.join(lines))

def addAtomSysName(project, ccpCode, atomName, altSysNames=(), sysName=None, 
                   atomSubType=1, molType='protein', namingSystem='CCPN_REMED'):
  """Add AtomSysName object
  """
  
  if sysName is None:
    sysName = atomName
  
  chemComp = genIo.getChemComp(project, molType, ccpCode)
  if chemComp is None:
    print '### ChemComp not found:', molType, ccpCode
  
  else:
    systemObj = chemComp.findFirstNamingSystem(name=namingSystem)
    if systemObj is None:
      systemObj = chemComp.newNamingSystem(name=namingSystem)
    
    nameObj = systemObj.findFirstAtomSysName(atomName=atomName, 
                                             atomSubType=atomSubType)
    
    if nameObj is None:
      systemObj.newAtomSysName(atomName=atomName, atomSubType=atomSubType,
                               sysName=sysName, altSysNames=altSysNames)
    else:
      print '### AtomSysName already exists:', molType, ccpCode, namingSystem, atomName, atomSubType
      
  

def createLogFile(identifier, task):
  
  dataDir = os.path.join(allDataDir, identifier[1:3], identifier)
  logDir = os.path.join(dataDir, 'log_' + task)
  if not os.path.isdir(logDir):
    os.makedirs(logDir)
  ff = '%s.%s.log' % (identifier,time.strftime("%Y.%m.%d_%H.%M"))
  logFileHandle = open(os.path.join(logDir, ff),'w')
  
  sys.stdout = logFileHandle
  
  return logFileHandle
  
  
def replacePdbNames(path):
  #def generalPdbFix(path):
  """ Fix PDB files, removing some errors. No effect on correct files
  """
  print '### fixing PDB', path
  
  replaces = (
  ('  CD  ILE',
   '  CD1 ILE'),
  (' 2CH  TRP',
   ' CH2  TRP'),
  (' 1NH  ARG',
   ' NH1  ARG'),
  (' 2NH  ARG',
   ' NH2  ARG'),
  ('  1HE2GLN',
   ' 1HE2 GLN'),
  ('  2HE2GLN',
   ' 2HE2 GLN'),
  ('  1HG1VAL',
   ' 1HG1 VAL'),
  ('  2HG1VAL',
   ' 2HG1 VAL'),
  ('  3HG1VAL',
   ' 3HG1 VAL'),
  ('  1HG2VAL',
   ' 1HG2 VAL'),
  ('  2HG2VAL',
   ' 2HG2 VAL'),
  ('  3HG2VAL',
   ' 3HG2 VAL'),
  ('  1HG1ILE',
   ' 1HG1 ILE'),
  ('  2HG1ILE',
   ' 2HG1 ILE'),
  ('  1HG2ILE',
   ' 1HG2 ILE'),
  ('  2HG2ILE',
   ' 2HG2 ILE'),
  ('  3HG2ILE',
   ' 3HG2 ILE'),
  ('  1HD1ILE',
   ' 1HD1 ILE'),
  ('  2HD1ILE',
   ' 2HD1 ILE'),
  ('  3HD1ILE',
   ' 3HD1 ILE'),
  ('  1HD1LEU',
   ' 1HD1 LEU'),
  ('  2HD1LEU',
   ' 2HD1 LEU'),
  ('  3HD1LEU',
   ' 3HD1 LEU'),
  ('  1HD2LEU',
   ' 1HD2 LEU'),
  ('  2HD2LEU',
   ' 2HD2 LEU'),
  ('  3HD2LEU',
   ' 3HD2 LEU'),
  ('  1HH1ARG',
   ' 1HH1 ARG'),
  ('  2HH1ARG',
   ' 2HH1 ARG'),
  ('  1HH2ARG',
   ' 1HH2 ARG'),
  ('  2HH2ARG',
   ' 2HH2 ARG'),
  ('  1HD2ASN',
   ' 1HD2 ASN'),
  ('  2HD2ASN',
   ' 2HD2 ASN'),
  ('  1HG2THR ',
   ' 1HG2 THR '),
  ('  2HG2THR ',
   ' 2HG2 THR '),
  ('  3HG2THR ',
   ' 3HG2 THR '),
   # NBNB the following two transform to MOLMOL names. 
   # But since these seem to wrok Ok for other TRP side chain protons
   # let us keep them
  (' 1HD  TRP',
   '  HD  TRP'),
  (' 2HH  TRP',
   '  HH  TRP'),
  )
  
  text = open(path).read()
  
  for src, dst in replaces:
    text = text.replace(src, dst)
  #
  open(path,'w').write(text)
  

def makeNmrCalcRun(ccpnProject, task='CCPN'):
  """ make NmrCalc.run in NmrCalcStore named 'task', 
  and add all possible data from NmrProject
  """
  
  taskRun = prepareNmrCalcRun(ccpnProject, task)
  
  # MolSystem 
  molSystem = ccpnProject.findFirstMolSystem()
  if molSystem:
    taskRun.newMolSystemData(name='firstMolSystem', 
                             details='auto added for %s' % task,
                             molSystem=molSystem)
    
    # Structures
    ll = molSystem.sortedStructureEnsembles()
    if ll:
      taskRun.newStructureEnsembleData(name='recentStructureEnsemble', 
                                       details='auto added for %s' % task,
                                       structureEnsemble=ll[-1])
    
  # NMR data
  nmrProject = taskRun.nmrCalcStore.nmrProject
  
  # Nmr Constaints
  ll = nmrProject.sortedNmrConstraintStores()
  if ll:
    taskRun.newConstraintStoreData(name='allRecentConstraints', 
                                   details='auto added for %s' % task,
                                   nmrConstraintStore=ll[-1])
  
  for nmrexp in nmrProject.sortedExperiments():
    
    # Spectra
    ll = nmrexp.sortedDataSources()
    if ll:
      spectrum = ll[-1]
      taskRun.newSpectrumData(name='allSpectra',
                              details='auto added for %s' % task,
                              dataSource=spectrum)
      
      # Peak Lists
      for peakList in spectrum.sortedPeakLists():
        if peakList.peaks:
          # Skip empty lists.
          taskRun.newPeakListData(name=peakList.name 
                                  or '%s:%s' % (spectrum.name, peakList.serial),
                                  details='auto added for %s' % task,
                                  peakList=peakList)
    
  # Shift and other measurement lists
  for ml in nmrProject.sortedMeasurementLists():
    if ml.measurements:
      # Skip empty lists.
      taskRun.newMeasurementListData(name=ml.name 
                                     or '%s:%s' % (ml.className, ml.serial),
                                     details='auto added for %s' % task,
                                     measurementList=ml)
      
      
  
def prepareNmrCalcRun(ccpnProject, task):
  """ set up ccpnProject to have an NmrCalc.run for the task
  """
  
  # get NmrCalcStore
  taskStore = ccpnProject.findFirstNmrCalcStore(name=task)
  if taskStore is None:
    nmrProject = (ccpnProject.currentNmrProject or 
                  ccpnProject.findFirstNmrProject() or 
                  ccpnProject.newNmrProject(name='auto'))
    taskStore = ccpnProject.newNmrCalcStore(name=task, nmrProject=nmrProject)
  
  taskRun = taskStore.newRun(status='pending', 
                             details='Autoprepared for %s' % task)
  #
  return taskRun


#######################################################################
#
# Gereral Util:
#
# NBNB should be moved elseqhere
#

'''
# Incomplete. May or may not be needed
def getFileInfo(fileDir, fileName):
  """estimate file type info for reading dta files 
  """
  result = {}
  
  head, ext = os.path.splitext(fileName.lower())
  
  result['name'] = head
  
  if ext == 'upl':
    result['dataType'] = 'distanceConstraints'
  
  elif ext == 'aco':
    result['dataType'] = 'dihedralConstraints'
  
  elif ext == 'tbl':
    
    if 'dihe' in head:
      result['dataType'] = 'dihedralConstraints'
    else:
      for tag in ('ambig', 'unambig', 'noe', 'dist', 'assigned', 'satisfied'):
        if tag in head:
          result['dataType'] = 'distanceConstraints'
          break
   '''       
      
  
  
    
    
  

  
if __name__ == '__main__':
  
  pass
  
  # Make raw json file from downloaded htm file of uploaded results
  #htmlFile = '/home/rhf22/Testing/CASD-NMR/results/unifi_it_targets.html'
  #extractToJson(htmlFile)
  
  #fixFourLetterHNames('/home/rhf22/Testing/CASD-NMR/results/patches/OR135/281/structures/or135_rosetta_ensemble.pdb')
  #for path in (
  # '/home/rhf22/Testing/CASD-NMR/results/patches/HR2876B/279/structures/HR2876B_r3.pdb',
  # '/home/rhf22/Testing/CASD-NMR/results/patches/HR6470A/135/structures/HR6470A_rdc.pdb',
  # '/home/rhf22/Testing/CASD-NMR/results/patches/HR6470A/136/structures/HR6470A.pdb',
  # '/home/rhf22/Testing/CASD-NMR/results/patches/YR313A/277/structures/YR313A_r3.pdb',
  # '/home/rhf22/Testing/CASD-NMR/results/patches/OR135/207/structures/OR135_cns.pdb',
  # '/home/rhf22/Testing/CASD-NMR/results/OR135/312/structures/OR135_cns.pdb', 
  # '/home/rhf22/Testing/CASD-NMR/results/HR2876B/313/structures/HR2876B_cns.pdb', 
  # '/home/rhf22/Testing/CASD-NMR/results/HR2876B/314/structures/HR2876B_r3.pdb', 
  # '/home/rhf22/Testing/CASD-NMR/results/YR313A/315/structures/YR313A_cns.pdb', 
  # '/home/rhf22/Testing/CASD-NMR/results/YR313A/316/structures/YR313A_r3.pdb', 
  
  #):
  #  replaceInSecondLines(path)
  #for path in (
  #   '/home/rhf22/Testing/CASD-NMR/results/patches/HR6470A/186/structures/final_auotNOE_hr6470_tails.pdb',
  #  ):
  #  replacePdbNames(path)
  
  
  
  
  
  
