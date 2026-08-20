from ccp.format.general import Util as ccpGenUtil
from memops.universal import Io as uniIo
import re, os

# next three obsolete?
from ccp.general.Constants import code1LetterToCcpCodeDict
from ccp.format.bruker.generalIO import BrukerGenericFile
from ccp.format.bruker.Constants import stdExpNameDict, pulProgFile

pulProgNameMap = None

specialCases = []
def parsePulProgName(name, fileLocation=None):
  """ Input: pulse program name to be parsed
             location of pulse program map file
             - if not given use standard file
      Output: List of (expName,isReversed) tuples
              expName can be either an NmrExpPrototype name 
              or a RefExperiment name.
              isReversed can be True, False, or None (meaning not set).
      
      Caches pulse program parsing map
  """
  
  global pulProgNameMap
  
  # get pulse program name translation map
  if fileLocation is None:
    if pulProgNameMap is None:
      fileLocation = uniIo.joinPath(uniIo.getTopDirectory(), 'data', 'ccp', 
                                    'bruker', 'BrukerPulseProgMap.xml')
      parseInfo = ccpGenUtil.loadAllParseInfo(fileLocation)
      pulProgNameMap = ccpGenUtil.convertParseInfo(parseInfo)
  
  else:
    parseInfo = ccpGenUtil.loadAllParseInfo(fileLocation)
    pulProgNameMap = ccpGenUtil.convertParseInfo(parseInfo)
  
  # get match result
  matchMap = matchName(name, pulProgNameMap, specialCases)
  dd = matchMap
  
  return matchMap['matches']
    

def testNameParsing(dirPath):
  """ Frame code to run name testing
  """
  
  # set of imput info
  stdInfoFile = uniIo.joinPath(uniIo.getTopDirectory(),
                               'data', 'ccp', 'bruker', 'BrukerPulseProgMap.xml')
  fileNames = os.listdir(dirPath)
  fileNames.sort()
  parseInfo = ccpGenUtil.loadAllParseInfo(stdInfoFile)
  useInfo = ccpGenUtil.convertParseInfo(parseInfo)
  
  # set up mapId:list dictionary
  matchByMap = {}
  notmatched = []
  skipped = []
  for expNameMap in parseInfo['expNameMaps']:
    matchByMap[id(expNameMap)] = []
  
  # match filenames and store by matchDict
  for fname in fileNames:
    dd = matchName(fname, useInfo, specialCases)
    dd['fname'] = fname
    expNameMap = dd.get('expNameMap')
    if (dd.get('status') == 'skip' 
        or not os.path.isfile(uniIo.joinPath(dirPath,fname))):
      skipped.append(dd)
    elif expNameMap is None:
      notmatched.append(dd)
    else:
      matchByMap[id(expNameMap)].append(dd)
    #print fname, dd['matches'], dd['names'], dd['unused']
  
  # print result
  print '\nSKIPPED:'
  for dd in skipped:
    print '    - ', dd['fname'], dd.get('usenames'), dd.get('prefixes'), dd.get('unused') 

  print '\nNO MATCH:'
  for dd in notmatched:
    print '    - ', dd['fname'], dd.get('usenames'), dd.get('prefixes'), dd.get('unused') 

  unused = []
  for xmap in parseInfo['expNameMaps']:
    ll = matchByMap[id(xmap)]
    prototypes = [(x.get('name'), x.get('isReversed')) for x in xmap['expPrototypes']]
    if ll:
      print '\nMATCH:', xmap.get('namesStartWith'), xmap.get('hasnames'), xmap.get('hasprefixes')
      print '      ', prototypes
      for dd in ll:
        if dd.get('unused'):
          ss = ' NB - '
        else:
          ss = '    - '
          
        print ss, dd['fname'], dd.get('usenames'), dd.get('prefixes'), dd.get('unused')
          
  for xmap in unused:
    print '\nUNUSED', xmap.get('namesStartWith'), xmap.get('hasnames'), xmap.get('hasprefixes') 
    print '      ', [(x.get('name'), x.get('isReversed')) for x in xmap['expPrototypes']]

  
def matchName(fname, useInfo, specialCases):
  
  result = ccpGenUtil.parseProgramName(fname, useInfo, specialCases)
  
  result['expNameMap'] = None
  result['matches'] = matches = []
  
  if result['names']:
    nameMaps = useInfo['nameMaps']
    names = result['names'] + result['endnames']
    names = result['usenames'] = [nameMaps.get(name, name) for name in names]
    nameset = set(names)
    prefixset = set(result['prefixes'])
    for dd in useInfo['expNameMaps']:
      starters = dd['namesStartWith']
      if names[:len(starters)] == starters:
        if nameset.issuperset(dd['hasnames']):
          if prefixset.issuperset(dd['hasprefixes']):
            
            for dd2 in dd['expPrototypes']:
              name = dd2.get('name')
              isReversed = dd2.get('isReversed')
              if isReversed in ('false', 'False'):
                isReversed = False
              elif isReversed in ('true', 'True'):
                isReversed = True
              matches.append((name, isReversed))
            
            result['expNameMap'] = dd
            break
  
  else:
    result['usenames'] =  result['endnames']
  
  #
  return result
  
  

def parseMusicExp(text, useInfo):
  """ parse music experiment of the form music_xyz_something
  treated as no prefixes, names == ['music', 'xyz'],
  parse 'something' for endnames
  """
  
  ll = text.split('_')
  
  result = {
   'names': ll[:2],
   'prefixes': [],
  }
  
  unused, endnames, suffix = ccpGenUtil.parseEndNames(useInfo['endnames'], 
                                                      '_'.join(ll[2:]))
  result['unused'] = unused
  result['endnames'] = endnames
  result['suffix'] = suffix
  if unused:
    result['status'] = 'partial'
  else:
    result ['status'] = 'OK'
  #
  return result

specialCases.append((re.compile("music_[a-zA-Z]+_"), parseMusicExp))

def parseRdExp(text, useInfo):
  """
  parse reduced-dimensionality expname of the form ".*rd_[a-zA-Z0-9]+_ij.*
  where ij are single-digit integers
  """
  
  indexChars = '123456789'
  ll = text.split('_')
  for ii in range(len(ll)-2):
    if ll[ii].endswith('rd'):
      dimsName = ll[ii+2][:2]
      if dimsName[0] in indexChars and dimsName[1] in indexChars:
        # we have a match
        if len(ll[ii+2]) == 2:
          del ll[ii+2]
        else:
          ll[ii+2] = ll[ii+2][2:]
          
        if len(ll[ii]) == 2:
          del ll[ii]
        else:
          ll[ii] = ll[ii][:-2]
        
        newText = '_'.join(ll)
        result = ccpGenUtil.parseProgramName(newText, useInfo)
        result['prefixes'].append('rd_')
        result['endnames'].insert(0,dimsName)
        #
        return result

specialCases.append((re.compile(".*rd_[a-zA-Z0-9]+_.*"), parseRdExp))

def getStdExpNamesByLength():

  """
  OBSOLETE?
  This gets the information about 'standard' Bruker experiment names when part
  of the pulse program name, and orders them by length.
  """

  stdExpNames = stdExpNameDict.keys()
  
  byLen = {}
  for stdExpName in stdExpNames:
    strLen = len(stdExpName)
    if not byLen.has_key(strLen):
      byLen[strLen] = []
    byLen[strLen].append(stdExpName)
    
  strLengths = byLen.keys()
  strLengths.sort()
  strLengths.reverse()
  
  stdExpNamesByLen = []
  
  for strLen in strLengths:
    byLen[strLen].sort()
    stdExpNamesByLen.extend(byLen[strLen])
  
  return stdExpNamesByLen
    
class PulProgInfo(BrukerGenericFile):

  """
  OBSOLETE?
  This reads the Bruker Pulprog.info file and extracts the Bruker codes that define
  aspects of a particular experiment (from the experiment name string).
  """

  def getInfo(self):
    
    #
    # These are added as default...
    #
    
    brukerCodeDict = {'2d': '2D experiment', '3d': '3D experiment', '4d': '4D experiment'}
    
    #
    # Start reading...
    #
    
    ppin = open(self.name)
    lines = ppin.readlines()
    ppin.close()
    
    for line in lines:

      if self.patt['emptyline'].search(line) or self.patt[self.format + 'SemiColonComment'].search(line):
        continue
      
      line = line.strip()
      
      brukerCode = None
      for spacing in ("     ","    "):
        if line.count(spacing):
          (brukerCode,info) = line.split(spacing)
          break

      if brukerCode:
        brukerCodeDict[brukerCode] = info

    return brukerCodeDict

  
class ParseBrukerExpNames:
  """OBSOLETE?
  """
  
  def __init__(self):
  
    self.stdExpNameDict = stdExpNameDict
    self.stdExpNamesByLen = getStdExpNamesByLength()
    
    ppi = PulProgInfo(pulProgFile)
    self.brukerCodeDict= ppi.getInfo()
    self.brukerCodes = self.brukerCodeDict.keys()
    self.brukerCodes.sort()
    
    for brukerCode in self.brukerCodes[:]:
      if len(brukerCode) == 3:
        self.brukerCodes.pop(self.brukerCodes.index(brukerCode))
        self.brukerCodes.insert(0,brukerCode)
    

  def loopBrukerExpNames(self,fileName):
  
    fin = open(fileName)
    lines = fin.readlines()
    fin.close()
    
    for line in lines:
    
      fileName = line.strip()
      
      if fileName[-5:] == '.info':
        continue
      elif fileName in ('README',):
        continue
       
      self.disectBrukerExpName(fileName)

  def disectBrukerExpName(self,brukerExpName):
  
    getAminoAcids = False

    print "%s:" % brukerExpName
    print

    #
    # First get main components, by largest string first
    #
    
    hasMainComponent = False
    for stdExpName in self.stdExpNamesByLen:

      if stdExpName[0] == '^':
        funcName = 'startswith'
        stdExpNameSearch = stdExpName[1:]
      else:
        funcName = 'count'
        stdExpNameSearch = stdExpName

      if getattr(brukerExpName,funcName)(stdExpNameSearch):
        print "  Main component: %s (%s)" % (stdExpNameSearch,self.stdExpNameDict[stdExpName])
        brukerExpName = brukerExpName.replace(stdExpNameSearch,'*' * len(stdExpNameSearch))
        hasMainComponent = True
        
        if stdExpName in ('music',):
          getAminoAcids = True
    
    if not hasMainComponent:
      print "  MAIN COMPONENT MISSING!"
          
    #
    # Get amino acids, if necessary. Assume these are between _, first _xxx_ block is amino acid one.
    #
    
    if getAminoAcids:
      nameComps = brukerExpName.split("_")
      aaNameComp = nameComps[1].upper()

      for code3Letter in code1LetterToCcpCodeDict['protein'].values() + ['Trpe']:
        if aaNameComp.count(code3Letter.upper()):
          print "    Bruker info: targeted at %s amino acid." % code3Letter
          aaNameComp = aaNameComp.replace(code3Letter.upper(),'*' * len(code3Letter))

      for code1Letter in code1LetterToCcpCodeDict['protein'].keys():
        if aaNameComp.count(code1Letter.upper()):
          print "    Bruker info: targeted at %s amino acid." % code1LetterToCcpCodeDict['protein'][code1Letter]
          aaNameComp = aaNameComp.replace(code1Letter.upper(),'*')

      brukerExpName = brukerExpName.replace(nameComps[1],aaNameComp.lower())

    #
    # Remove _
    #
    
    while brukerExpName.count("_"):
      brukerExpName = brukerExpName.replace("_","*")

    #
    # Then get Bruker components (two to three letter code). Do this from back of string.
    #

    # TODO: READ FROM BACK, FIRST THROW AWAY NUMBERS (BUT CHECK IF NOT ONE OF THE MEANINGFUL ONES)
    # Try first by throwing numbers away, see if works, then try using numbers?
    #
    # Otherwise do some special casing here - rl1, rl2 end on ineptrl1, ineptrl2
    
    for strInd in range(len(brukerExpName) - 1,1,-1):
      if brukerExpName[strInd] != '*' and brukerExpName[strInd - 1] != '*':
        for brukerCode in self.brukerCodes:
          if brukerExpName[strInd-len(brukerCode)+1:strInd+1].count(brukerCode):
            rIndex = brukerExpName.rfind(brukerCode)
            brukerExpName = brukerExpName[:rIndex] + brukerExpName[rIndex:].replace(brukerCode,'*' * len(brukerCode))
            print "    Bruker info: %s (%s)" % (brukerCode,self.brukerCodeDict[brukerCode])
    
    #
    # Remove .1, .2 type stuff
    #
    
    for i in range(1,10):
      for prefix in ('.',''):
        endCode = '%s%d' % (prefix,i)
        if brukerExpName.endswith(endCode):
          brukerExpName = brukerExpName.replace(endCode, '*' * len(endCode))
          
    #
    # Remove numbers left in name...
    #
    
    for i in range(1,10):
      if brukerExpName.count("*%d*" % i):
        brukerExpName = brukerExpName.replace("*%d*" % i,"***")
     
    #
    # If 'i' - means intra. Check if left
    #
   
    isIntra = False 
    if brukerExpName.startswith("i*"):
      brukerExpName = brukerExpName.replace('i*','**')
      isIntra = True
    elif brukerExpName.count("*i*"):
      brukerExpName = brukerExpName.replace('*i*','***')
      isIntra = True
   
    if isIntra:
      print "    Bruker info: intra experiment."
            
    #
    # Print out name if anything not recognized...
    #
    
    if len(brukerExpName) != brukerExpName.count('*'):
      print
      print "  Remaining unparsed name: %s" % brukerExpName

    print 
