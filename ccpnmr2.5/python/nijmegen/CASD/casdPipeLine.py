""" Commands for CASD analysis pipeline
"""


import os, sys, json, re
import tarfile, zipfile, bz2, gzip, shutil, tempfile
import traceback

from nijmegen.CASD import Constants as casdConstants
from nijmegen.CASD import Util as casdUtil
from memops.general import Io as genIo

from ccp.lib import StructureIo

from pdbe.deposition.dataFileImport.formatConverterWrapper import FormatConverterWrapper

#casdNmrDir = casdConstants.casdNmrDir
allDataDir = casdConstants.allDataDir

#casdResultsDir = casdConstants.casdResultsDir
#casdInputCcpnDir = casdConstants.casdInputCcpnDir
#casdInputDir = casdConstants.casdInputDir
#resultPatchesDir = os.path.join(casdResultsDir, 'patches')

import codecs, locale
sys.stdout = codecs.getwriter(locale.getpreferredencoding())(sys.__stdout__) 

def restraintOverview(entryNames, extractDir=None):
    """ Make overview  of restraint data available, and check data type
    """
   
    result = {}
   
    if not extractDir:
      extractDir = tempfile.mkdtemp(dir=casdConstants.topTmpDir)
    
    for entryName in entryNames:
      result[entryName] = info = {}
      inputDir = os.path.join(allDataDir, entryName[1:3], entryName, 
                        entryName+'.input', 'restraints')
      ll = os.listdir(inputDir)
      if len(ll) == 1:
        source = os.path.join(inputDir, ll[0])
      
        # restraints, check them
        targetDir = os.path.join(extractDir, entryName[1:3], entryName)
        if not os.path.exists(targetDir):
          casdUtil.extractCompressedFile(source, targetDir, entryName)
        checkDir = casdUtil.getLowestSubDir(targetDir, 
                                            followDirs=('cns_format',))
        
        convertType = checkRestraintTypes(checkDir, entryName)
      
      elif ll:
        print entryName, 'ERROR, multifiles', ll
      
      else:
        print entryName, 'NONE'

def checkRestraintTypes(restraintDir, entryName=None):
  """ count restraints of differnet types and return restraint reader type
  """
  extentionTypes = {
   '.upl':'CYANA',
   '.lol':'CYANA',
   '.aco':'CYANA',
   '.str':'NMR-STAR',
   '.tbl':'CNS',
  }
  
  
  counter = {}
  
  for fname in os.listdir(restraintDir):
    ff,fext = os.path.splitext(fname)
    ftype = extentionTypes.get(fext, '???')
    ll = counter.get(ftype)
    if ll:
      ll.append(fname)
    else:
      counter[ftype] = [fname]
  
  ll =  [(9999, entryName or os.path.basename(restraintDir), str(len(counter)))]
  for tag, val in sorted(counter.items()):
    ll.append((len(val), tag, val))
  
  ss = ''
  for tt in reversed(sorted(ll)):
    ss += "%s %s; " % (tt[1], tt[2])
  
  print ss 


def makeOverview(resultData, fieldOrder):
  """ make list-of-lists of results data, sorted.
  """
  
  
  result = []
  for info in resultData:
    
    entryName = casdUtil.getEntryName(info)
    
    # get normal data
    row = []
    result.append(row)
    for tag in fieldOrder:
      if tag == 'EntryName':
        row.append(entryName)
      else:
        row.append(info.get(tag))
    
    # Get data-is-present code string
    ll = []
    pp = os.path.join(allDataDir, entryName[1:3], entryName, entryName)
    if not os.path.exists(pp+'.tgz'):
      ll.append('Ccpn_NO')
    
    inputDir = pp + '.input'
    if not os.listdir(os.path.join(inputDir, 'restraints')):
      ll.append('Restraints_NO')
    if not os.listdir(os.path.join(inputDir, 'structures')):
      ll.append('Structures_NO')
    if os.listdir(os.path.join(inputDir, 'superseded')):
      ll.append('UseOrigData_NO')
    row.append(' '.join(ll))
    
  #
  result.sort()
  return  result

def makeOverview2(calcData, fieldOrder):
  """ make list-of-lists of results data, sorted.
  """
  
  
  result = []
  for entryName, info in calcData.items():
    
    # get normal data
    row = []
    result.append(row)
    for tag in fieldOrder:
      if tag == 'EntryName':
        row.append(entryName)
      else:
        row.append(info.get(tag))
    
    # Get data-is-present code string
    ll = []
    pp = os.path.join(allDataDir, entryName[1:3], entryName, entryName)
    if not os.path.exists(pp+'.tgz'):
      ll.append('Ccpn_NO')
    
    inputDir = pp + '.input'
    if not os.listdir(os.path.join(inputDir, 'restraints')):
      ll.append('Restraints_NO')
    if not os.listdir(os.path.join(inputDir, 'structures')):
      ll.append('Structures_NO')
    if os.listdir(os.path.join(inputDir, 'superseded')):
      ll.append('UseOrigData_NO')
    row.append(' '.join(ll))
    
  #
  result.sort()
  return  result

def countCategories(resultData, tags):
  """ take in a list of dictionaries, and classify by the values efined by tags 
  """
  result = {}
  for dd in resultData:
    tt = tuple(dd.get(tag) for tag in tags)
    count = result.get(tt,0)
    result[tt] = count + 1
  #
  return result

def forAllEntries(calcData, func, entries=()):
  """Select data dict and execute func
  """
  
  if entries:
    dicts = [calcData.get(x) for x in entries if x in calcData]
  else:
    dicts = [x for x in calcData.values() if not x.get('isOriginal')]
  
  print 'Executing %s %s times' % (func.__name__, len(dicts))
  for dd in dicts:
    func(dd)
  

def makeCalcData(inputData, resultData):
  """make calcData data structure from merging inputData and resultData
  """
  origNames = {}
  calcData = {}
  
  # get data from input 
  for dd in inputData:
    if not dd.get("Invalid"):
      dd2 ={'isOriginal':True}
      for tag in ('Target', 'LACS CA/CB Offset', 'Oligomeric state', 'PDBcode',
                   'Defined Residues', ):
        val = dd.get(tag)
        if val:
          val = str(val)
        dd2[tag] = val
      entryName = casdUtil.getEntryName(dd2, isOriginal=True)
      calcData[entryName] = dd2
      origNames[dd['Target']]= entryName
  
  # get data from entries
  for dd in resultData:
    if not dd.get("Invalid"):
      origName = origNames.get(dd['Target'])
      if origName:
        dd2 = {}
        
        dd2['EntryID'] = dd.get('EntryID')
        for tag in ('Target', 'Group', 'Program Type', 
                     "RDCdata", "Peaklist", "Truncated", 'Submitted on'):
          val = dd.get(tag)
          if val:
            val = str(val)
          dd2[tag] = val
          
        origDd = calcData[origName]
        for tag in ('LACS CA/CB Offset', 'Oligomeric state', 'PDBcode',
                    'Defined Residues', ):
          dd2[tag] = origDd[tag]
          
        entryName = casdUtil.getEntryName(dd2)
        calcData[entryName] = dd2
  #
  return calcData

def makeCcpnProject(entryName):
  """ Execute conversion to CCPN project
  """
  
  logFileHandle = None
  
  try:
  
    #entryName = casdUtil.getEntryName(info)
    orgName = entryName.split('_')[0] + '_Org'
    
    print 'Starting', entryName
    
    logFileHandle = casdUtil.createLogFile(entryName, 'extractEntry')
  
    # get CCPN project from Org data
    #orgName = casdUtil.getEntryName(info, isOriginal=True)
    path = os.path.join(allDataDir, orgName[1:3], orgName)
    ppath = os.path.join(path, orgName)
    ff = ppath + '.tgz'
    if not (os.path.exists(ppath) or os.path.exists(ff)):
      raise Exception("NO original CCPN project in %s" % path)
    if not os.path.isdir(ppath):
      casdUtil.extractCompressedFile(ff, path, entryName)
    if not os.path.exists(ppath):
      raise Exception("NO extracted CCPN project in %s" % path)
    ccpnProject = genIo.loadProject(ppath, suppressGeneralDataDir=True)
    
    
    # neutralize any other pending CASD-NMR projects
    casdRun = casdUtil.prepareNmrCalcRun(ccpnProject, 'CING')
    for run in casdRun.nmrCalcStore.findAllRuns(status='pending'):
      if run is not casdRun:
        run.status = 'provisional'
 
    #
    dataDir = os.path.join(allDataDir, entryName[1:3], entryName)
    tmpdir = tempfile.mkdtemp(dir=casdConstants.topTmpDir)
    try:
 
      # Extract structure data
      tmpstruc = os.path.join(tmpdir, 'structures')
      src = casdUtil.getInputFile(entryName, 'structures')
      casdUtil.extractCompressedFile(src, tmpstruc, entryName, okExts=('pdb',))
      tmpstruc = casdUtil.getLowestSubDir(tmpstruc, followDirs=('cns_format',))
      structureFiles = os.listdir(tmpstruc)
 
 
      # Extract restraint data
      tmprestr = os.path.join(tmpdir, 'restraints')
      src = casdUtil.getInputFile(entryName, 'restraints', ignoreErrors=True)
      if src:
        casdUtil.extractCompressedFile(src, tmprestr, entryName)
        tmprestr = casdUtil.getLowestSubDir(tmprestr, followDirs=('cns_format',))
        restraintFiles = os.listdir(tmprestr)
      else:
        restraintFiles = ()
        print 'WARNING, %s no restraints found at %s' % (entryName, src)
 
      # read in data
      
      # FormatConverter version
      fcw = FormatConverterWrapper(ccpnProject=ccpnProject)
      # dataIo version
      #fcw = None
      
      if structureFiles:
        # NBNB uses Rasmus in-development trunk structure reading.
        # WOrks well. Temporarily disabled
        # read in structures
        pdbFiles = [x for x in structureFiles
                    if any(x.endswith(y) for y in casdConstants.pdbEndings)]
 
        floatFiles = [x for x in structureFiles if x.endswith('.float')]
 
        if floatFiles:
          # Use only pdb files with names that match float files
          stems = set(x[:-6] for x in floatFiles)
          pdbFiles = [x for x in pdbFiles if x[:-4] in stems]
 
        pdbPaths = [os.path.join(tmpstruc,x) for x in pdbFiles]
        
        if True:
        #if fcw is None:
	      #Always use dataIo version
          # dataIo version
          ensemble = StructureIo.getStructureFromFiles(
                                 ccpnProject.findFirstMolSystem(), pdbPaths)
                                             
          if ensemble is None:
            print '### Skipping %s, no structures loaded' % entryName
            
          else:
            print '### num files, ensemble', len(pdbPaths), ensemble.ensembleId
            casdRun.newStructureEnsembleData(name=entryName, 
                                             structureEnsemble=ensemble)
        
        else:
          # FormatConverter version
          #fileInfo = fcw.determineFileInfo(pdbPaths[0])
          if len(pdbPaths) != 1:
            print 'WARNING %s pdb files, only one read. TBD FIX' % len(pdbPaths)
          dataType = 'coordinates'
          formatName = 'pseudoPdb'
          pdbPath = pdbPaths[0]
          print 'Reading structure file', dataType, formatName, pdbPath
          fcw.readFile(dataType, formatName, pdbPath)
            
          
          # NBNB TODO 1) How to set up trying true PDB before pseudoPdb?
          #           2) How to read several files into an ensemble?
          #           3) How to get hold of the new ensemble for putting in NmrCalc
        
      else:
        print '### Skipping %s, no structure file' % entryName
        
    
      # Make NmrCalc object for shift list
      # NBNB consider later: if we are reading in assigned peaks, 
      # shifts may change. NBNB
      shiftLists = casdRun.nmrCalcStore.nmrProject.findAllMeasurementLists(className='ShiftList')
      if len(shiftLists) == 1:
        casdRun.newMeasurementListData(name='Shiftlist', 
                                       measurementList=shiftLists.pop())
      else:
        print 'WARNING. %s shift lists found, should be s' % len(shiftLists)
    
      # Restraints reading 
      if restraintFiles:
        if fcw is None:
          # NBNB TBD dataIo restraint reading to go here
          #for rfile in restraintFiles:
          #  fileInfo = casdUtil.getFileInfo(tmprestr, rfile)
          pass
        
        else:
          # FormatConverter version
          restraintLists = []
          for rfile in restraintFiles:
            rpath = os.path.join(tmprestr,rfile)
            fileInfo = fcw.determineFileInfo(rpath)
            dataType = fileInfo.get('dataType')
            formatName = fileInfo.get('formatName')
            if dataType is None or formatName is None:
              print 'Skipping unidentified restraint file', dataType, formatName,  rfile
            
            elif dataType not in ('distanceConstraints', 'dihedralConstraints',
                                  'rdcConstraints',):
              print 'Skipping wrong type of restraint file', dataType, formatName,  rfile
              
            else:
              print 'Reading restraint file', dataType, formatName, rfile
              fcw.readFile(dataType, formatName, rpath)
              if fcw.conversionSuccess:
                print ("Successful restraint file read:\n%s" % fcw.conversionInfo)
                restraintLists.append(fcw.ccpnObjectOrList)
              else:
                print ("Failed restraint file read:\n%s" % fcw.conversionInfo)

          if restraintLists:
            print ("Found restraint lists: %s" % len(restraintLists))
            casdRun.newConstraintStoreData(constraintLists=restraintLists, name='Restraintlists')

        
        # linkResonances
        print '### linking resonances'
        linkingInfo = fcw.linkAllResonancesToAtoms()

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
        pass
    
    # rename and package project
    ccpnOutputDir = os.path.join(dataDir, entryName)
    genIo.saveProject(ccpnProject, newPath=ccpnOutputDir, 
              newProjectName=entryName, checkValid=True, removeExisting=True)
    genIo.packageProject(ccpnProject, ccpnOutputDir)
    shutil.rmtree(ccpnOutputDir)
    ccpnOutputPath = ccpnOutputDir + '.tgz'
    print 'SUCCESS, %s saved to %s' % (entryName, ccpnOutputPath)
    
    return ccpnOutputPath
  
  except:
    print 'ERROR for %s' % (entryName)
    traceback.print_exc(file=sys.stdout)
  
  finally:
    if logFileHandle is not None:
      logFileHandle.close() 
      sys.stdout = sys.__stdout__

'''
def extractResults(infoDict):
  """ Extract result data 
  """
  
  # set up
  target = infoDict['Target']
  entry = infoDict['EntryID']
  dataId = '%s %s' % (target, entry)
  print '### ', dataId
  
  targetDir = os.path.join(casdNmrDir, 'results', target, str(entry))
  structDir = os.path.join(targetDir, 'structures')
  structFile = os.path.join(targetDir, infoDict['Structure'])
  restraintDir = os.path.join(targetDir, 'restraints')
  
  # extract structure files
  casdUtil.extractCompressedFile(structFile, structDir, dataId, okExts=('.pdb',))
  
  # get subdir containing data to load
  loadDir = casdUtil.getLowestSubDir(structDir)
  #loadDir = casdUtil.getLowestSubDir(structDir, followDirs=('cns_format',))
  
  # check for missing / unexpected structure files
  ll = [x for x in os.listdir(loadDir)]
  ll2 = [x for x in ll if not x.endswith('.pdb')]
  if ll2:
    print ('### WARNING, %s non-pdb structure data: %s' %
           (dataId, ll2))
  if not ll or len(ll2) == len(ll):
    print ('### WARNING, %s NO pdb structure data' % dataId)
    
  # Make always-applicable PDB fixes
  #for (eachDir, junk, files) in os.walk(loadDir):
  #  ll = [x for x in files if x.endswith('.pdb')]
  #  for ff in ll:
  #    generalPdbFix(os.path.join(eachDir, ff))
  
  # copy patches across to structures dir
  # NBNB now disabled. Files are fixed up front instead.
  #patchesDir = os.path.join(resultPatchesDir, target, str(entry), 'structures')
  #if os.path.isdir(patchesDir):
  #  ll = os.listdir(patchesDir)
  #  for ff in ll:
  #    print '### structure patch:', dataId, ff
  #    shutil.copy(os.path.join(patchesDir, ff), loadDir)
  
  # extract restraint files
  if infoDict['Restraints']:
    restraintFile = os.path.join(targetDir, infoDict['Restraints'])
    casdUtil.extractCompressedFile(restraintFile, restraintDir, dataId + ' restraints')

  # copy patches across to restraints dir
  # NBNB now disabled. Files are fixed up front instead.
  #loadDir = casdUtil.getLowestSubDir(structDir, followDirs=('cns_format',))
  #patchesDir = os.path.join(resultPatchesDir, target, str(entry), 'restraints')
  #if os.path.isdir(patchesDir):
  #  ll = os.listdir(patchesDir)
  #  for ff in ll:
  #    print '### restraints patch:', dataId, ff
  #    shutil.copy(os.path.join(patchesDir, ff), loadDir)
  '''

'''
# Superseded version
def makeCcpnProject(info):
  """ Execute conversion to CCPN project
  """
  
  ccpnOutputPath = None
    
  try:
    target = info['Target']
    entry = info['EntryID']
    projId = '%s_%s' % (target,entry)
    
    if entry in casdConstants.skipEntries:
      print '### Skipping %s, marked to ignore' % projId
      return
 
    # Temporary - must be improved later
 
    ccpnInputDir = os.path.join(casdNmrDir, 'input', 'ccpn', target)
    targetDir = os.path.join(casdNmrDir, 'results', target, str(entry))
    structureDir = os.path.join(targetDir, 'structures')
    #structureDir = casdUtil.getLowestSubDir(structureDir, followDirs=('cns_format',))
    structureDir = casdUtil.getLowestSubDir(structureDir)
    restraintDir = os.path.join(targetDir, 'restraints')
    restraintDir = casdUtil.getLowestSubDir(restraintDir)
    #restraintDir = casdUtil.getLowestSubDir(restraintDir, followDirs=('cns_format',))
    ccpnOutputPath = os.path.join(targetDir, 'ccpn')
 
    restraintFiles = os.listdir(restraintDir)
    structureFiles = os.listdir(structureDir)
 
    if structureFiles:
      # Structure only
      pdbFiles = [x for x in structureFiles
                  if any(x.endswith(y) for y in casdConstants.pdbEndings)]
      
      floatFiles = [x for x in structureFiles if x.endswith('.float')]
      
      if floatFiles:
        # Use only pdb files with names that match float files
        stems = set(x[:-6] for x in floatFiles)
        pdbFiles = [x for x in pdbFiles if x[:-4] in stems]
 
      ccpnProject = genIo.loadProject(ccpnInputDir,
                                      suppressGeneralDataDir=True)
 
      pdbPaths = [os.path.join(structureDir,x) for x in pdbFiles]
      StructureIo.getStructureFromFiles(ccpnProject.findFirstMolSystem(), pdbPaths)
      
      print '### num files, ensembles', len(pdbPaths), len(ccpnProject.structureEnsembles)
      x = ccpnProject.findFirstStructureEnsemble()
      print '### num structures: %s' % (x and len(x.models))
 
    else:
      print '### Skipping %s, no structure files' % projId
      return
 
    # rename and package project
    if not os.path.exists(ccpnOutputPath):
      os.makedirs(ccpnOutputPath)
    ccpnOutputPath = os.path.join(ccpnOutputPath, projId)
    genIo.saveProject(ccpnProject, newPath=ccpnOutputPath, 
              newProjectName=projId, checkValid=True, removeExisting=True)
    genIo.packageProject(ccpnProject, ccpnOutputPath)
    shutil.rmtree(ccpnOutputPath)
    ccpnOutputPath += '.tgz'
    #genIo.saveProject(ccpnProject, ccpnOutputPath, newProjectName=projId,
    #                  removeExisting=True, checkValid=True,
    #                  changeDataLocations=True)
    print 'SUCCESS, %s %s saved to %s' % (projId, info['Program'], 
                                          ccpnOutputPath)
    
    return ccpnOutputPath
    
  except:
    print 'ERROR for %s %s' % (projId, info['Program'])
    traceback.print_exc(file=sys.stdout)
    '''
'''
def expandResultData(resultData):
  """ Hack to add pdb codes, project names, and Program Types to resultData.json
  """
  # Program type dict:
  progTypes = {
   ('AnnArbor','a'):'I-TASSER',
   ('Cheshire','a'):'Cheshire',
   ('Cheshire','b'):'Cheshire-YAPP',
   ('Frankfurt','a'):'CYANA',
   ('Lyon','a'):'UNIO',
   ('Madison','a'):'Ponderosa',
   ('Munich','a'):'CS-DP-Rosetta',
   ('Paris','a'):'ARIA',
   ('Piscataway','a'):'ASDP-CNS',
   ('Piscataway','b'):'ASDP-Rosetta',
   ('Seattle','a'):'CS-HM-Rosetta',
   ('Seattle','b'):'CS-HM-DP-Rosetta',
   ('Trieste','a'):'BE-metadynamics	',
   ('Utrecht','a'):'CS-Rosetta',
   ('Utrecht','b'):'CS-DP-Rosetta',
  }
  
  
  # make PDB code dict
  jsonFileIn = os.path.join(casdNmrDir, 'input', 'inputData.json')
  inputData = json.load(open(jsonFileIn))
  pdbCodes = {}
  for dd in inputData:
    target = dd['Target (download)']
    pdbCodes[target] = dd['PDB ID'].lower()
  
  for dd in resultData:
    entryID = dd['EntryID']
    if entryID > 120:
      group = dd['Group']
      dd['Program Type'] = progTypes.get((group,dd['Subgroup']))
      pdbCode = pdbCodes.get(dd['Target'])
      dd['PDBcode'] = pdbCode
      dd['EntryName'] = '%s_%s_%s' % (pdbCode, group, entryID)
'''

'''
def oldToNewData(info):
  """ Hack - make new data tree and copy data across
  """
    
  # make new directories
  entryName = casdUtil.getEntryName(info)
  target = info['Target']
  topdir = entryName[1:3]
  
  inputDir = os.path.join(casdConstants.allDataDir, topdir, entryName, 
                          entryName+'.input')
  if not os.path.isdir(inputDir):
    os.makedirs(inputDir)
  
  newdirs = {}
  for ss in ('structures', 'restraints', 'superseded'):
    newdirs[ss] = path = os.path.join(inputDir, ss)
    if not os.path.isdir(path):
      os.mkdir(path)
  
  if entryName.endswith('Org'):
    # This is starting data
    
    # get restraints
    useDir = os.path.join(casdNmrDir, 'input', 'downloads')
    ll = [x for x in os.listdir(useDir) if x.startswith(target)]
    if len(ll) == 1:
      shutil.copy(os.path.join(useDir,ll[0]), newdirs['restraints'])
    else:
      print 'ERROR, found %s restraint files for ' % len(ll), entryName
    
    # get pdb file
    pdbCode = info['PDBcode'].lower()
    useDir = os.path.join(casdNmrDir, 'input', 'pdb')
    ll = [x for x in os.listdir(useDir) if x and x.lower().startswith(pdbCode)]
    if len(ll) == 1:
      shutil.copy(os.path.join(useDir,ll[0]), newdirs['structures'])
    else:
      print 'ERROR, found %s restraint files for ' % len(l), entryName
  
  
  else:
    # this is a results sumbmission
    oldDataDir = os.path.join(casdNmrDir, 'results', info['Target'], str(info['EntryID']))
    
    for junk1, junk2, files in os.walk(oldDataDir):
      # just a quick way of getting files-only contents
      break
    
    # Structure file
    ff = info.get('Structure')
    if ff in files:
      shutil.copy(os.path.join(oldDataDir,ff), newdirs['structures'])
      files.remove(ff)
    
    # Restraints file
    ff = info.get('Restraints')
    if ff in files:
      shutil.copy(os.path.join(oldDataDir,ff), newdirs['restraints'])
      files.remove(ff)
    
    # Superseded files
    for ff in files:
      shutil.copy(os.path.join(oldDataDir,ff), newdirs['superseded'])
  
  print '### Copied', entryName, target
  '''


if __name__ == '__main__':
  
  inputFile = os.path.join(allDataDir, 'inputData.json')
  inputData = json.load(open(inputFile))
  
  resultFile = os.path.join(allDataDir, 'resultData.json')
  resultData = json.load(open(resultFile))
  
  calcDataFile = os.path.join(allDataDir, 'calcData.json')
  calcData = makeCalcData(inputData, resultData)
  json.dump(calcData, open(calcDataFile, 'w'), sort_keys=True, indent=4)
  
  '''
  #expandResultData(resultData)
  #json.dump(resultData, open(resultFile, 'w'), sort_keys=True, indent=4)
  '''
  #tags = ('Target', 'PDBcode')
  #tags = ('Program Type','Group','Peaklist','PDBcode','RDCdata', 'Truncated')
  #tags = ('RDCdata','Program Type',)
  #dd = countCategories(calcData.values(), tags)
  #print tags
  #for tag,val in sorted(dd.items()):
  #  print val, tag
  # extract results to CCPN projects
  #entries = []
  #doAdd = False    
  #for dd in resultData:
  #  indx = dd['EntryID']
  #  if doAdd:
  #    entries.append(indx)
  #  elif indx == 291:
  #    doAdd = True
  
  #restraintOverview(sorted(calcData.keys()), '/home/rhf22/tmpdata/tmp2T11BX')
  
  #entryList = ['2ltm_Piscataway_279']
  #entryList = ['2m5o_Paris_303']
  #entryList = ['2m5o_Piscataway_301']
  entryList = sorted(list(x for x in calcData.keys() if not x.endswith('_Org')))
  #entryList = ['2la6_Paris_144', '2l9r_Piscataway_135', '2m5o_Paris_303']

  for entryName in entryList:
    makeCcpnProject(entryName)
  
  #for entryName in ('2lah_Cheshire_321','2lah_Cheshire_322','2lah_Cheshire_323',
  #                  '2lci_Cheshire_324','2ltm_Cheshire_325','2ltm_Cheshire_326',
  #                  '2ltl_Cheshire_327','2ltl_Cheshire_328','2m2e_Cheshire_329',
  #):
  #for entryName in [x for x in calcData if '2m5o_' in x and not 'Org' in x]:
  #   makeCcpnProject(entryName)
  
  # Make result overview file
  fieldOrder = ('Group', 'Program Type', 'EntryName', 'Peaklist', 'RDCdata', 
                'Truncated', "Submitted on", 'Target')
  #overview = makeOverview2(calcData, fieldOrder)
  #overview = makeOverview((x for x in resultData 
  #                         if x['EntryID'] > 120 and not x.get('Invalid')),
  #                        fieldOrder)
  #print fieldOrder
  #for ll in overview:
  #  print '\t'.join(unicode(x) for x in ll)
  
  # counting classes:
  #for cat,count in sorted(countCategories(resultData, ('Program Type',)).items()):
  #  print count, ':'.join(unicode(x) for x in cat)
  
  
  # Directory analysis:
  #count = 0
  #for info in resultData:
  #  if info['EntryID'] > 120 and not info.get('Invalid'):
  #    count += 1
  #    entryName = casdUtil.getEntryName(info)
  #    pp = os.path.join(allDataDir, entryName[1:3], entryName)
  #    ll = os.listdir(pp)
  #    ll2 = os.listdir(os.path.join(pp, 'log_extractEntry'))
  #    if (len(ll) != 3 or not (entryName + '.input' in ll and 
  #                             entryName + '.tgz' in ll and 
  #                             'log_extractEntry' in ll)
  #        or not ll2):
  #      print 'WARNING, incorrect files:', entryName, ll
  #print 'ENTRIES: ', count
