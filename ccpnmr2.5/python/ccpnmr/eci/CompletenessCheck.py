import re, sys

from ccpnmr.eci.AditMandFields import aditMandDict

from memops.general.Io import loadProject

from pdbe.nmrStar.IO import nmrStarDict
from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport
from pdbe.nmrStar.IO.Ccpn_To_NmrStar import Ccpn_To_NmrStar

# TODO: entry keywords, struc genomics, constraints, software vendors/tasks

CCPN_OBJ_NAME_DICT_SHORT = {'nmrEntry':                'BMRB Entry',
                            'contactPersons':          'Contact Person(s)',
                            'personInGroups':          'Contact Person Address',
                            'authors':                 'Deposition Authors',
                            #'sgInfo':                  'Structural Genomics',
                            #'keywords':                'Structural Keywords',
                            'primaryCitation':         'Main Citation',
                            'otherCitations':          'Other Citations',
                            'citationAuthors':         'Citation Authors',
                            'molSystem':               'Molecular System',
                            'chains':                  'Molecular System Chains',
                            'molecule':                'Macromolecule',
                            'entryMolecules':          'Biological Sources',
                            'naturalSource':           'Natural Source',
                            'experimentalSource':      'Experimental Source',
                            'sample':                  'Sample',
                            'sampleComponent':         'Sample Components',
                            'sampleConditionSet':      'Sample Condition Set',
                            'sampleConditions':        'Sample Conditions',
                            'software':                'Software',
                            #'vendorName':              'Software Vendor',
                            #'tasks':                   'Software Tasks',
                            'experiments':             'NMR Experiments',
                            'spectrometer':            'NMR Spectrometer',
                            'spectrometerManu':        'Spectrometer Manufacturer',
                            'shiftReferences':         'Chemical Shift References',
                            'shiftLists':              'Chemical Shift Lists',
                            'shiftExperiments':        'Shift List Experiments',
                            'shiftSoftware':           'Shift List Software',
                            'shifts':                  'Chemical Shifts',
                            'structureGenerations':    'Structural Calculation Data',
                            'strucGenSoftware':        'Structure Generation Software',
                            'structureEnsemble':       'Structural Ensemble',
                            'models':                  'Ensemble Models',
                            'nmrConstraintStore':      'NMR Constraint Store',
                            'distanceConstraintLists': 'Distance Constraint Lists',
                            'distanceConstraints':     'Distance Constraints',}

# TODO: add more detail to this dictionary

CCPN_OBJ_NAME_DICT_LONG = {'nmrEntry':                'BMRB Entry',
                           'contactPersons':          'Contact Person(s)',
                           'personInGroups':          'Contact Person Address',
                           'authors':                 'Deposition Authors',
                           #'sgInfo':                  'Structural Genomics',
                           #'keywords':                'Structural Keywords',
                           'primaryCitation':         'Main Citation',
                           'otherCitations':          'Other Citations',
                           'citationAuthors':         'Citation Authors',
                           'molSystem':               'Molecular System',
                           'chains':                  'Molecular System Chains',
                           'molecule':                'Macromolecule',
                           'entryMolecules':          'Biological Sources',
                           'naturalSource':           'Natural Source',
                           'experimentalSource':      'Experimental Source',
                           'sample':                  'Sample',
                           'sampleComponent':         'Sample Components',
                           'sampleConditionSet':      'Sample Condition Set',
                           'sampleConditions':        'Sample Conditions',
                           'software':                'Software',
                           #'vendorName':              'Software Vendor',
                           #'tasks':                   'Software Tasks',
                           'experiments':             'NMR Experiments',
                           'spectrometer':            'NMR Spectrometer',
                           'spectrometerManu':        'Spectrometer Manufacturer',
                           'shiftReferences':         'Chemical Shift References',
                           'shiftLists':              'Chemical Shift Lists',
                           'shiftExperiments':        'Shift List Experiments',
                           'shiftSoftware':           'Shift List Software',
                           'shifts':                  'Chemical Shifts',
                           'structureGenerations':    'Structural Calculation Data', 
                           'strucGenSoftware':        'Structure Generation Software',
                           'structureEnsemble':       'Structural Ensemble',
                           'models':                  'Ensemble Models',
                           'nmrConstraintStore':      'NMR Constraint Store',
                           'distanceConstraintLists': 'Distance Constraint Lists',
                           'distanceConstraints':     'Distance Constraints',}


STAR_TO_CCPN_SF_NAMES = {'assembly':              'Molecular System',
                         'citations':             'Citations',
                         'entry_interview':       'BMRB Entry Data',
                         'entity':                'Macromolecule',
                         'sample':                'Sample',
                         'sample_conditions':     'Sample Condition Set',
                         'experiment_list':       'NMR Experiment List',
                         'NMR_spectrometer_list': 'NMR Spectrometer List',
                         'experimental_source':   'Experimental Source List',
                         'entry_information':     'BMRB Entry',
                         'natural_source':        'Natural Source List',
                         'chem_comp':             'Chemical Components',
                         'software':              'Software'}

STAR_TO_CCPN_TABLE_NAMES = {'Entity_assembly':           'Molecular System Chains',
                            'Citation_author':           'Citation Authors',
                            'Sample_component':          'Sample Components',
                            'Sample_condition_variable': 'Sample Conditions',
                            'Experiment':                'NMR Experiments',
                            'NMR_spectrometer_view':     'NMR Spectrometer',
                            'Entity_experimental_src':   'Experimental Source',
                            'Contact_person':            'Contact Person(s)',
                            'Entry_author':              'Deposition Author',
                            'Entity_natural_src':        'Natural Source',
                            'Task':                      'Software Task',
                            'Vendor':                    'Software Vendor/Author'}


def getAttrOrFunc(obj,attrOrFunc):

  patt = re.compile("(\(.*\))")

  returnObj = None

  if obj:

    funcSearch = patt.search(attrOrFunc)

    if funcSearch:
      funcArgs = funcSearch.group(1)
      funcName = attrOrFunc.replace(funcArgs,'')

      func = getattr(obj,funcName)
      returnObj = eval("func%s" % funcArgs)

    else:
      returnObj = getattr(obj,attrOrFunc)

  return returnObj

def getCcpnObject(ccpnMap, startObject):

  ccpnObjectStrings = ccpnMap.split('.')

  findObject = None

  if len(ccpnObjectStrings) > 1:
    try:
      findObject = reduce(lambda obj, attrOrFunc: getAttrOrFunc(obj,attrOrFunc),ccpnObjectStrings[1:],startObject)

    except:
      print "  Cannot get CCPN object for link %s, starting object %s" % (ccpnObjectStrings, startObject)
      raise
      findObject = []

  else:
    findObject = startObject

  return findObject

def getCcpnMapValue(ccpnMap, startObject):

  ccpnMapValue = None

  if type(ccpnMap) == type(''):
    ccpnMapValue = getCcpnObject(ccpnMap, startObject)

  elif type(ccpnMap) == type( () ):
    (ccpnVarName,getFunc) = ccpnMap

    try:
      ccpnMapValue = getFunc(getCcpnObject(ccpnVarName, startObject) )

    except:
      print "  Cannot CCPN map value for %s, function %s" % (ccpnVarName,getFunc)
      raise

  elif ccpnMap == None:
    pass

  else:
    try:
      ccpnMapValue = ccpnMap()

    except:
      self.setMessage("  Error: ccpnLoop can only be string, function or tuple(string,function). Is %s." % str(ccpnMap) )

  return ccpnMapValue

def getSfCcpnMapValues(ccpnMap, sfKey, tagName):

  if type(ccpnMap) == type(''):
    ccpnMapString = ccpnMap

  elif type(ccpnMap) == type( () ):
    (ccpnMapString, getFunc) = ccpnMap

  ccpnObjectStrings = ccpnMapString.split('.')

  topMap = ccpnObjectStrings[0]
  endMap = '.'.join(ccpnObjectStrings)

  if not sfKey in sfCcpnMapDict:
    sfCcpnMapDict[sfKey] = {}

  if not tagName in sfCcpnMapDict[sfKey]:
    sfCcpnMapDict[sfKey][tagName] = set()

  if topMap == 'nmrEntry':
    if tagName not in ('Dep_release_code_nmr_exptl', 'Dep_release_code_sequence'):
      startObject = ccpnVar[topMap]
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      sfCcpnMapDict[sfKey][tagName].add(ccpnValue)
    else:
      del(sfCcpnMapDict[sfKey][tagName])

  elif topMap == 'citation':
    startObject = ccpnVar['primaryCitation']

    citationClass = getCcpnMapValue('citation.className', startObject)

    journalCount = 0

    if tagName != 'Journal_abbrev':
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      sfCcpnMapDict[sfKey][tagName].add(ccpnValue)
    elif citationClass == 'JournalCitation' and tagName == 'Journal_abbrev':
      journalCount += 1
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      sfCcpnMapDict[sfKey][tagName].add(ccpnValue)

    for startObject in ccpnVar['otherCitations']:
      
      citationClass = getCcpnMapValue('citation.className', startObject)

      if tagName != 'Journal_abbrev':
        ccpnValue = getCcpnMapValue(ccpnMap, startObject)
        sfCcpnMapDict[sfKey][tagName].add(ccpnValue)
      elif citationClass == 'JournalCitation' and tagName == 'Journal_abbrev':
        journalCount += 1
        ccpnValue = getCcpnMapValue(ccpnMap, startObject)
        sfCcpnMapDict[sfKey][tagName].add(ccpnValue)

    # TODO: test this code somehow.

    startObject = ccpnVar['primaryCitation']

    if not journalCount and startObject and not sfCcpnMapDict[sfKey][tagName]:
      del(sfCcpnMapDict[sfKey][tagName])

  elif topMap == 'molecule':
    for chain in ccpnVar['chains']:
      startObject = chain.molecule
      if tagName == 'Polymer_seq_one_letter_code':
        #molType = getCcpnObject('molecule.molType', startObject)
        seqLength = getCcpnObject('molecule.seqLength', startObject)
        isStdLinear = getCcpnObject('molecule.isStdLinear', startObject)
        if seqLength == 1 or not isStdLinear:
          continue

      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      sfCcpnMapDict[sfKey][tagName].add(ccpnValue)

  elif topMap == 'software':
    for startObject in ccpnVar['software']:
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      sfCcpnMapDict[sfKey][tagName].add(ccpnValue)

  elif topMap == 'sample':
    for expt in ccpnVar['experiments']:
      startObject = expt.sample
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      sfCcpnMapDict[sfKey][tagName].add(ccpnValue)

  # TODO:
  # - No tags in Eldon's ADIT-NMR dictionary for shift lists.

  elif topMap in ('shiftList',):

    #print 'SF: [%s] [%s]' % (topMap, ccpnMap)
    return False

  else:
    startObject = ccpnVar[topMap]
    ccpnValue = getCcpnMapValue(ccpnMap, startObject)
    sfCcpnMapDict[sfKey][tagName].add(ccpnValue)
    #print 'ELSE SF: [%s] [%s] [%s]' % (topMap, ccpnMap, ccpnValue)

  return True

def getTableCcpnMapValues(ccpnMap, sfKey, tableName, tagName):

  if type(ccpnMap) == type(''):
    ccpnMapString = ccpnMap

  elif type(ccpnMap) == type( () ):
    (ccpnMapString, getFunc) = ccpnMap

  ccpnObjectStrings = ccpnMapString.split('.')

  topMap = ccpnObjectStrings[0]
  endMap = '.'.join(ccpnObjectStrings)

  #print 'TOP: [%s] END: [%s]' % (topMap, endMap)

  if not sfKey in tableCcpnMapDict:
    tableCcpnMapDict[sfKey] = {}

  if not tableName in tableCcpnMapDict[sfKey]:
    tableCcpnMapDict[sfKey][tableName] = {}

  if not tagName in tableCcpnMapDict[sfKey][tableName]:
    tableCcpnMapDict[sfKey][tableName][tagName] = set()

  if topMap == 'contactPerson':
    for startObject in ccpnVar['contactPersons']:
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

    if 'personInGroups' in ccpnObjDict and not ccpnObjDict['personInGroups'][0]:
      #print 'LIST OF TAGS: [%s]' % starDict['entry_information']['tables']['Contact_person']['tagNames'] # Minus name attrs.
      if tagName in ('Email_address', 'Address_1', 'City', 'State_province', 'Country', 'Postal_code', 'Phone_number', 'Organization_type', 'Role'):
        del(tableCcpnMapDict[sfKey][tableName][tagName])

  elif topMap == 'author':
    for startObject in ccpnVar['authors']:
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  elif topMap == 'chain':
    for startObject in ccpnVar['chains']:
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  elif topMap == 'entryMolecule' and tagName != 'Vector_type':
    for startObject in ccpnVar['entryMolecules']:
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  elif topMap == 'naturalSource':
    for entryMol in ccpnVar['entryMolecules']:
      startObject = entryMol.molecule.naturalSource
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  elif topMap == 'experimentalSource' or (topMap == 'entryMolecule' and tagName == 'Vector_type'):
    exptSourceCount = 0

    for entryMol in ccpnVar['entryMolecules']:

      prodMethod = getCcpnMapValue('entryMolecule.productionMethod', entryMol)

      if prodMethod != 'recombinant technology':
        if tagName in ('Host_org_scientific_name', 'Vector_type', 'Vector_name'):
          continue
      elif tagName == 'Host_org_scientific_name':
        exptSourceCount += 1
          
      startObject = entryMol.experimentalSource
      if tagName == 'Vector_type':
        startObject = entryMol

      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

      # TODO: vectorType still needs to be set in the CCPN project.

      #print 'DATA: [%s] [%s] [%s] [%s] [%s] [%s] [%s]' % (sfKey, tableName, tagName, startObject, ccpnMap, ccpnValue, prodMethod)

    # TODO: this needs testing

    if not exptSourceCount and not tableCcpnMapDict[sfKey][tableName][tagName]:
      del(tableCcpnMapDict[sfKey][tableName][tagName])

  elif topMap == 'software':
    for startObject in ccpnVar['software']:
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue) 

  elif topMap == 'sampleComponent':
    for expt in ccpnVar['experiments']:
      sample = expt.sample
      for startObject in sample.sampleComponents:
        ccpnValue = getCcpnMapValue(ccpnMap, startObject)
        tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  elif topMap == 'sampleCondition':
    for expt in ccpnVar['experiments']:
      scs = expt.sampleConditionSet
      for startObject in scs.sampleConditions:
        ccpnValue = getCcpnMapValue(ccpnMap, startObject)
        tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  # TODO: this will be linked to nmrEntry soon.

  elif topMap == 'spectrometer':
    for expt in ccpnVar['experiments']:
      startObject = expt.spectrometer
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  elif topMap == 'experiment':
    for startObject in ccpnVar['experiments']:
      ccpnValue = getCcpnMapValue(ccpnMap, startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  elif topMap == 'task':
    for startObject in ccpnVar['software']:
      ccpnValue = getCcpnMapValue('software.tasks', startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(', '.join(ccpnValue) )

  elif topMap == 'vendorName':
    for startObject in ccpnVar['software']:
      ccpnValue = getCcpnMapValue('software.vendorName', startObject)
      tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)

  # TODO:
  # - Keywords (and structural genomics) only required for a structural deposition.
  # - No tags in Eldon's ADIT-NMR dictionary for shift references.
  # - Too complicated to do measurementByIndividualAtom?

  elif topMap in ('keyword', 'chemShiftRef', 'measurementByIndividualAtom'):

    #print 'TABLE: [%s] [%s]' % (topMap, ccpnMap)
    return False

  else:
    startObject = ccpnVar[topMap]
    ccpnValue = getCcpnMapValue(ccpnMap, startObject)
    tableCcpnMapDict[sfKey][tableName][tagName].add(ccpnValue)
    #print 'ELSE TABLE: [%s] [%s] [%s]' % (topMap, ccpnMap, ccpnValue)

  return True

def setCcpnVar(nmrEntry):

  global ccpnVar

  ccpnVar = {'nmrEntry': nmrEntry}

  for ccpnLink in ccpnLinkList:
    ccpnMap = 'nmrEntry.%s' % ccpnLink
    ccpnVar[ccpnLink] = getCcpnObject(ccpnMap, nmrEntry)

  ccpnVar['chains'] = getCcpnObject('molSystem.chains', nmrEntry.molSystem)

  #print 'CCPNVAR: [%s]' % ccpnVar

def checkSfTags(sfKey):

  sfTagDict = {}

  if sfKey not in starDict or 'tagNames' not in starDict[sfKey]:
    return sfTagDict

  for tagName in starDict[sfKey]['tagNames']:
    #mandTagFlag = starDict[sfKey]['tags'][tagName][3]
    mandTagFlag = True

    if mandTagFlag:
      if tagName in ('ID', 'Sf_category', 'Entry_ID'):
        continue

      if sfKey not in sfTagDict:
        sfTagDict[sfKey] = []

      sfTagDict[sfKey].append(tagName)

  #print 'SF TAG DICT: [%s] [%s]' % (sfKey, sfTagDict)

  return sfTagDict

def checkMappingSfTags(sfTagDict):

  sfMappingDict = {}

  for sfKey in sfTagDict.keys():
    if sfKey not in sfMappingDict:
      sfMappingDict[sfKey] = {}

    for tagName in sfTagDict[sfKey]:

      if tagName in ccpn2StarDict[sfKey]['tags']:
        tagMapping = ccpn2StarDict[sfKey]['tags'][tagName]

        if tagName not in sfMappingDict[sfKey]:
          ccpnMapValue = getSfCcpnMapValues(tagMapping, sfKey, tagName)

          if ccpnMapValue:
            sfMappingDict[sfKey][tagName] = tagMapping

          #print 'SF TAG: [%s] [%s] [%s] [%s]' % (sfKey, tagName, tagMapping, ccpnMapValue)

      elif 'CONDITIONAL' in ccpn2StarDict[sfKey]['tags']:
        for key1 in ccpn2StarDict[sfKey]['tags']['CONDITIONAL']:
          for key2 in ccpn2StarDict[sfKey]['tags']['CONDITIONAL'][key1]:
            if tagName in ccpn2StarDict[sfKey]['tags']['CONDITIONAL'][key1][key2]:
              tagMapping = ccpn2StarDict[sfKey]['tags']['CONDITIONAL'][key1][key2][tagName]

              if tagName not in sfMappingDict[sfKey]:
                ccpnMapValue = getSfCcpnMapValues(tagMapping, sfKey, tagName)

                if ccpnMapValue:
                  sfMappingDict[sfKey][tagName] = tagMapping

                #print 'COND SF TAG: [%s] [%s] [%s] [%s]' % (sfKey, tagName, tagMapping, ccpnMapValue)

      elif (tagName + '=LOCAL') in ccpn2StarDict[sfKey]['tags']:
        pass
        #print 'SF LOCAL TAG: [%s] [%s] [%s]' % (
        #  sfKey, tagName, ccpn2StarDict[sfKey]['tags'][tagName + '=LOCAL'])

      # Always set 'Sf_framecode' in NmrStarExport - so no mapping in ccpn2StarDict.

      elif tagName == 'Sf_framecode':
        pass

      else:
        pass # TODO: add me
        #print 'NO SF TAG: [%s] [%s]' % (sfKey, tagName)

  return sfMappingDict

def checkTableTags(sfKey, tableName):

  tableTagDict = {}

  for tagName in starDict[sfKey]['tables'][tableName]['tagNames']:
    #mandTagFlag = starDict[sfKey]['tables'][tableName]['tags'][tagName][3]
    mandTagFlag = True

    if mandTagFlag:
      if tagName in ('ID', 'Sf_category', 'Entry_ID'):
        continue

      if sfKey not in tableTagDict:
        tableTagDict[sfKey] = {}

      if tableName not in tableTagDict[sfKey]:
        tableTagDict[sfKey][tableName] = []

      tableTagDict[sfKey][tableName].append(tagName)

  #print 'TABLE TAG DICT: [%s] [%s] [%s]' % (sfKey, tableName, tableTagDict)

  return tableTagDict

def checkMappingTableTags(tableTagDict):

  tableMappingDict = {}

  for sfKey in tableTagDict.keys():
    if sfKey not in tableMappingDict:
      tableMappingDict[sfKey] = {}

    for tableName in tableTagDict[sfKey].keys():
      if tableName not in tableMappingDict[sfKey]:
        tableMappingDict[sfKey][tableName] = {}

      if tableName not in ccpn2StarDict[sfKey]['tables']:
        # TODO: add me
        #print 'NO TABLE: [%s] [%s]' % (sfKey, tableName)
        continue

      for tagName in tableTagDict[sfKey][tableName]:

        if tagName in ccpn2StarDict[sfKey]['tables'][tableName]['tags']:
          tagMapping = ccpn2StarDict[sfKey]['tables'][tableName]['tags'][tagName]
          if tagName not in tableMappingDict[sfKey][tableName]:
            ccpnMapValue = getTableCcpnMapValues(tagMapping, sfKey, tableName, tagName)

            if ccpnMapValue:
              tableMappingDict[sfKey][tableName][tagName] = tagMapping

            #print 'TABLE TAG: [%s] [%s] [%s] [%s] [%s]' % (sfKey, tagName, tableName, tagMapping, ccpnMapValue)

        elif (tagName + '=LOCAL') in ccpn2StarDict[sfKey]['tables'][tableName]['tags']:
          pass
          #print 'TABLE LOCAL TAG: [%s] [%s] [%s]' % (
          #  sfKey, tagName, ccpn2StarDict[sfKey]['tables'][tableName]['tags'][tagName + '=LOCAL'])

        else:
          pass # TODO: add me
          #print 'NO TABLE TAG: [%s] [%s] [%s]' % (sfKey, tableName, tagName)

  return tableMappingDict

def checkCcpnObjects(nmrEntry, submissionType):

  entryObj = 'nmrEntry'

  ccpnObjDict = {entryObj: (True, CCPN_OBJ_NAME_DICT_SHORT[entryObj], CCPN_OBJ_NAME_DICT_LONG[entryObj])}

  ccpnLinkList.append('primaryCitation')
  ccpnLinkList.append('otherCitations')
  #ccpnLinkList.append('laboratories')
  #ccpnLinkList.append('chains')
  ccpnLinkList.append('software')

  contAddrFlag = True

  #sgInfoFlag = True
  #keywordFlag = True

  citAuthFlag = True

  chainFlag = True
  molFlag = True

  natSourceFlag = True
  expSourceFlag = True

  #vendorFlag = True
  #taskFlag = True

  scsFlag = True
  scFlag = True

  sampleFlag = True
  sampCompFlag = True

  nmrSpecFlag = True
  nmrManuFlag = True

  shiftRefFlag = True

  shiftListFlag = True
  shiftFlag = True
  shiftExptFlag = True
  shiftSoftFlag = True

  strucGenSoftFlag = True

  strucEnsFlag = True
  strucModelsFlag = True

  constStoreFlag = True
  distConstListFlag = True
  distConstFlag = True

  for ccpnLink in ccpnLinkList:
    if ccpnLink in ('study',):
      continue

    ccpnObjs = getattr(nmrEntry, ccpnLink)

    #print 'CCPN OBJ: [%s] [%s]' % (ccpnLink, ccpnObjs)

    if not ccpnObjs:
      #print 'No CCPN objects for %s' % ccpnLink
      ccpnObjDict[ccpnLink] = (False, CCPN_OBJ_NAME_DICT_SHORT[ccpnLink], CCPN_OBJ_NAME_DICT_LONG[ccpnLink])

    else:
      ccpnObjDict[ccpnLink] = (True, CCPN_OBJ_NAME_DICT_SHORT[ccpnLink], CCPN_OBJ_NAME_DICT_LONG[ccpnLink])

    if ccpnLink == entryObj:
      pass

      #if not ccpnObjDict[ccpnLink][0]:
      #  sgInfoFlag = keywordFlag = False

      #if not ccpnObjs.sgInfo:
      #  sgInfoFlag = False

      #if not ccpnObjs.keywords:
      #  keywordFlag = False

    elif ccpnLink == 'contactPersons':
      if not ccpnObjDict[ccpnLink][0]:
        contAddrFlag = False

      for contactPerson in ccpnObjs:
        if not contactPerson.findFirstPersonInGroup():
          contAddrFlag = False

    elif ccpnLink == 'primaryCitation':
      if not ccpnObjs:
        citAuthFlag = False
        continue

      if not ccpnObjDict[ccpnLink][0]:
        citAuthFlag = False

      if not ccpnObjs.authors:
        citAuthFlag = False

    elif ccpnLink == 'otherCitations':
      for citation in ccpnObjs:
        if not citation.authors:
          citAuthFlag = False

    elif ccpnLink == 'molSystem':
      if not ccpnObjs:
        chainFlag = molFlag = False
        continue

      if not ccpnObjDict[ccpnLink][0]:
        chainFlag = molFlag = False

      if not ccpnObjs.chains:
        chainFlag = False

      for chain in ccpnObjs.chains:
        if not chain.molecule:
          molFlag = False
          
      # TODO: Could check for residues?

    elif ccpnLink == 'entryMolecules':
      if not ccpnObjDict[ccpnLink][0]:
        natSourceFlag = expSourceFlag = False

      for entryMolecule in ccpnObjs:
        if not entryMolecule.molecule.naturalSource:
          natSourceFlag = False
        if not entryMolecule.experimentalSource:
          expSourceFlag = False

    elif ccpnLink == 'software':
      pass

      #if not ccpnObjDict[ccpnLink][0]:
      #  vendorFlag = taskFlag = False

      #for software in ccpnObjs:
      #  if not software.vendorName:
      #    vendorFlag = False
      #  if not software.tasks:
      #    taskFlag = False

    elif ccpnLink == 'experiments':
      if not ccpnObjDict[ccpnLink][0]:
        scsFlag = scFlag = sampleFlag = sampCompFlag = nmrSpecFlag = nmrManuFlag = shiftRefFlag = False

      for experiment in ccpnObjs:
        if not experiment.sampleConditionSet:
          scsFlag = scFlag = False
        elif not experiment.sampleConditionSet.sampleConditions:
          scFlag = False

        if not experiment.sample:
          sampleFlag = sampCompFlag = False
        elif not experiment.sample.sampleComponents:
          sampCompFlag = False

        if not experiment.spectrometer:
          nmrSpecFlag = nmrManuFlag = False
        elif not experiment.spectrometer.manufacturer:
          nmrManuFlag = False

        if not experiment.shiftReferences and experiment.shiftList: # Only need reference if expt has a shift list
          shiftRefFlag = False

    elif ccpnLink == 'structureGenerations' and submissionType=='PDB':
      if not ccpnObjDict[ccpnLink][0]:
        strucGenSoftFlag = strucEnsFlag = strucModelsFlag = constStoreFlag = distConstListFlag = distConstFlag = False

      for structureGeneration in ccpnObjs:
        if not structureGeneration.method or not structureGeneration.method.software:
          strucGenSoftFlag = False
        if not structureGeneration.structureEnsemble:
          strucEnsFlag = strucModelsFlag = False
        elif not structureGeneration.structureEnsemble.models:
          strucModelsFlag = False

        if not structureGeneration.nmrConstraintStore:
          constStoreFlag = False
        else:
          distConstLists = structureGeneration.nmrConstraintStore.findAllConstraintLists(className='DistanceConstraintList')

          if not distConstLists:
            distConstListFlag = distConstFlag = False

          else:
            for distConstList in distConstLists:
              if not distConstList.constraints:
                distConstFlag = False

        # TODO: Could check for models, chains, residues, atoms, coords?

  shiftLists = nmrEntry.findAllMeasurementLists(className='ShiftList', isSimulated=False)

  shiftListFlag = shiftFlag = shiftExptFlag = shiftSoftFlag = True

  if not shiftLists:
    shiftListFlag = shiftFlag = shiftExptFlag = shiftSoftFlag = False

  else:
    for shiftList in shiftLists:
      if not shiftList.measurements:
        shiftFlag = False
      if not shiftList.experiments:
        shiftExptFlag = False
      if not shiftList.method or not shiftList.method.software:
        shiftSoftFlag = False

  objOtherNameDict = {'personInGroups':          contAddrFlag,
                      #'sgInfo':                  sgInfoFlag,
                      #'keywords':                keywordFlag,
                      'citationAuthors':         citAuthFlag,
                      'chains':                  chainFlag,
                      'molecule':                molFlag,
                      'naturalSource':           natSourceFlag,
                      'experimentalSource':      expSourceFlag,
                      'sampleConditionSet':      scsFlag,
                      'sampleConditions':        scFlag,
                      'sample':                  sampleFlag,
                      'sampleComponent':         sampCompFlag,
                      'spectrometer':            nmrSpecFlag,
                      'spectrometerManu':        nmrManuFlag,
                      #'vendorName':              vendorFlag,
                      #'tasks':                   taskFlag,
                      'shiftLists':              shiftListFlag,
                      'shiftExperiments':        shiftExptFlag,
                      'shiftSoftware':           shiftSoftFlag,
                      'shifts':                  shiftFlag,
                      'shiftReferences':         shiftRefFlag}

  if submissionType == 'PDB':
    pdbObjOtherNameDict = {'structureEnsemble':       strucEnsFlag,
                           'models':                  strucModelsFlag,
                           'strucGenSoftware':        strucGenSoftFlag,
                           'nmrConstraintStore':      constStoreFlag,
                           'distanceConstraintLists': distConstListFlag,
                           'distanceConstraints':     distConstFlag}

    objOtherNameDict.update(pdbObjOtherNameDict)

  for objName in objOtherNameDict.keys():
    ccpnObjDict[objName] = (objOtherNameDict[objName], CCPN_OBJ_NAME_DICT_SHORT[objName], CCPN_OBJ_NAME_DICT_LONG[objName])

  return ccpnObjDict

def checkNmrEntryCompleteness(nmrEntry, textReport=False, submissionType='PDB'):

  global starDict
  global ccpn2StarDict
  global ccpnLinkList

  global ccpnObjDict

#  starDict = nmrStarDict.sfDict

  starDict = aditMandDict

  starExport = NmrStarExport(nmrEntry)

  ccpn2Star = Ccpn_To_NmrStar(starExport)
  ccpn2StarDict = ccpn2Star.sfDict

  ccpnLinkList = ccpn2Star.nmrEntryLinkList

  # These are set in ccpnMap/ccpnLoop routines in Ccpn_To_NmrStar

  ccpn2Star.nmrEntry = nmrEntry
  ccpn2Star.authors = list(nmrEntry.authors)

  #print 'DICTS: [%s] [%s]' % (starDict, ccpn2StarDict)
 
  ccpnObjDict = checkCcpnObjects(nmrEntry, submissionType)

  setCcpnVar(nmrEntry)

  #print 'CCPN OBJ DICT: [%s]' % ccpnObjDict

  finalInfoDict = {}

  redColour    = '#FF4040'
  orangeColour = '#E0B040'
  yellowColour = '#F0F020'
  greenColour  = '#B0FFB0'

  for ccpnObj in ccpnObjDict.keys():
    #print 'OBJ: [%s] [%s]' % (ccpnObj, ccpnObjDict[ccpnObj])

    if ccpnObj == 'otherCitations':
      continue

    if ccpnObj == 'structureGenerations' and submissionType == 'BMRB':
      continue

    if not ccpnObjDict[ccpnObj][0]:
      #print 'MISSING OBJ: [%s]' % ccpnObj
      finalInfoDict[ccpnObj] = (ccpnObjDict[ccpnObj][1], 'No "' + ccpnObjDict[ccpnObj][2] + '" set/selected', redColour, False, 'Object', 0, ccpnObj, None, None)
    else:
      #print 'PRESENT OBJ: [%s]' % ccpnObj
      finalInfoDict[ccpnObj] = (ccpnObjDict[ccpnObj][1], '"' + ccpnObjDict[ccpnObj][2] + '" done', greenColour, True, 'Object', 3, ccpnObj, None, None)


  # TODO list:

  # CCPN objects not in NMR-STAR
  #'entryMolecules',

  # CCPN objects that are not mandatory for NMR depositions - but required for structure depositions
  #'structureEnsemble',
  #'structureGenerations',
  # Constraint objects

  # NMR-STAR SF's/tables to be done - These aren't objects in CCPN - but are attributes.
  #'SG_project', 'Struct_keywords',
  #'Vendor', 'Task',


  # Rough translation of NmrStar save frames/tables to CCPN objects.

  starSf2CcpnDict = {'entry_information':        'nmrEntry',
                     'citations':                'primaryCitation', # 'otherCitations',
                     'assembly':                 'molSystem',
                     'entity':                   'molecule',
                     'natural_source':           'naturalSource',
                     'experimental_source':      'experimentalSource',
                     'sample':                   'sample',
                     'sample_conditions':        'sampleConditionSet',
                     'software':                 'software',
                     'NMR_spectrometer_list':    'spectrometer',
                     'experiment_list':          'experiments',
                     'assigned_chemical_shifts': 'shiftLists',
                     'chem_shift_reference':     'shiftReferences'}

  starTable2CcpnDict = {'Contact_person':            'contactPersons', # 'personInGroups',
                        'Entry_author':              'authors',
                        #'SG_project':                None,
                        #'Struct_keywords':           'keywords',
                        'Citation_author':           'citationAuthors',
                        'Entity_assembly':           'chains',
                        'Entity_natural_src':        'naturalSource',
                        'Entity_experimental_src':   'experimentalSource',
                        'Sample_component':          'sampleComponent',
                        'Sample_condition_variable': 'sampleConditions',
                        #'Vendor':                    'software',
                        #'Task':                      'software',
                        'NMR_spectrometer_view':     'spectrometer',
                        'Experiment':                'experiments',
                        'Chem_shift_experiment':     'shiftExperiments',
                        'Atom_chem_shift':           'shifts',
                        'Chem_shift_ref':            'shiftReferences'}

  # Mandatory save frames/tables in NmrStar (ADIT-NMR).

  # TODO: Eldon's ADIT-NMR dictionary doesn't include chem shift references,
  #   chem shift lists and structural (PDB) data.

  mandSfList = ['entry_information',
                'citations',
                'assembly',
                'entity',
                'natural_source',
                'experimental_source',
                'sample',
                'sample_conditions',
                'software',
                'NMR_spectrometer_list',
                'experiment_list',
                'chem_shift_reference',
                'assigned_chemical_shifts']

  mandTableDict = {'entry_information':        ['Contact_person',
                                                'Entry_author',
                                                'SG_project',
                                                'Struct_keywords'],
                   'citations':                ['Citation_author'],
                   'assembly':                 ['Entity_assembly'],
                   'entity':                   [],
                   'natural_source':           ['Entity_natural_src'],
                   'experimental_source':      ['Entity_experimental_src'],
                   'sample':                   ['Sample_component'],
                   'sample_conditions':        ['Sample_condition_variable'],
                   'software':                 ['Vendor',
                                                'Task'],
                   'NMR_spectrometer_list':    ['NMR_spectrometer_view'],
                   'experiment_list':          ['Experiment'],
                   'chem_shift_reference':     ['Chem_shift_ref'],
                   'assigned_chemical_shifts': ['Chem_shift_experiment',
                                                'Atom_chem_shift']}

  # These two dictionaries store the CCPN data.

  global sfCcpnMapDict
  global tableCcpnMapDict

  sfCcpnMapDict = {}
  tableCcpnMapDict = {}

  sfTagDict = {}

  for sfKey in nmrStarDict.sfList:
    if sfKey not in mandSfList:
      continue

    if sfKey in starSf2CcpnDict:
      ccpnObj = starSf2CcpnDict[sfKey]

      #print 'SF KEY: [%s] [%s] [%s]' % (sfKey, ccpnObj, ccpnObjDict[ccpnObj])

      # If no CCPN object, then we don't need to check for missing attributes.
      # Need to tell people to make a new CCPN object first.

      if ccpnObj in ccpnObjDict and not ccpnObjDict[ccpnObj][0]:
        #print 'NO SF CCPN OBJ: [%s] [%s]' % (sfKey, ccpnObj)
        continue

    sfTagDict = checkSfTags(sfKey)
    sfMappingDict = checkMappingSfTags(sfTagDict)

    # TODO: add me
    #print 'MAPPING SF DICT: [%s] [%s]' % (sfKey, sfMappingDict)

    if sfKey not in starDict or not 'tableNames' in starDict[sfKey]:
      continue

    for tableName in starDict[sfKey]['tableNames']:
      if tableName not in mandTableDict[sfKey]:
        continue

      if tableName in starTable2CcpnDict:
        ccpnObj = starTable2CcpnDict[tableName]

        #print 'TABLE KEY: [%s] [%s] [%s] [%s]' % (sfKey, tableName, ccpnObj, ccpnObjDict[ccpnObj])

        # Again if no CCPN object...

        if ccpnObj in ccpnObjDict and not ccpnObjDict[ccpnObj][0]:
          if ccpnObj == 'contactPersons':
            if 'personInGroups' in ccpnObjDict and not ccpnObjDict['personInGroups'][0]:
              #print 'NO TABLE CCPN OBJ: [%s] [%s/personInGroups]' % (sfKey, tableName, ccpnObj)
              continue
          else:
            #print 'NO TABLE CCPN OBJ: [%s] [%s]' % (sfKey, tableName, ccpnObj)
            continue

      tableTagDict = checkTableTags(sfKey, tableName)
      tableMappingDict = checkMappingTableTags(tableTagDict)

      # TODO: add me
      #print 'MAPPING TABLE DICT: [%s] [%s] [%s]' % (sfKey, tableName, tableMappingDict)

  for sfKey in sfCcpnMapDict:
    for tagName in sfCcpnMapDict[sfKey]:
      ccpnData = sfCcpnMapDict[sfKey][tagName]
      starInfo = starDict[sfKey]['tags'][tagName][3].strip()

      #print 'SF DATA: [%s] [%s] [%s] [%s]' % (sfKey, tagName, ccpnData, starInfo)

      emptyFlag = True

      if None in ccpnData:
        ccpnData.remove(None)
      elif '' in ccpnData:
        ccpnData.remove('')
      else:
        emptyFlag = False

      if not ccpnData:
        finalInfoDict[sfKey + '_' + tagName] = (STAR_TO_CCPN_SF_NAMES[sfKey] + ' [' + tagName + ']', 'No "' + starInfo + '" attribute set/selected', orangeColour, False, 'SF Attr', 1, sfKey, None, tagName)
      elif emptyFlag:
        finalInfoDict[sfKey + '_' + tagName] = (STAR_TO_CCPN_SF_NAMES[sfKey] + ' [' + tagName + ']', 'Some "' + starInfo + '" attributes not set/selected', yellowColour, False, 'SF Attr', 2, sfKey, None, tagName)
      #else:
      #  finalInfoDict[sfKey + '_' + tagName] = (sfKey + ': ' + tagName, starInfo + ' attribute set/selected', greenColour, True, 'SF Attr')

  for sfKey in tableCcpnMapDict:
    for tableName in tableCcpnMapDict[sfKey]:
      for tagName in tableCcpnMapDict[sfKey][tableName]:
        ccpnData = tableCcpnMapDict[sfKey][tableName][tagName]
        starInfo = starDict[sfKey]['tables'][tableName]['tags'][tagName][3].strip()

        #print 'TABLE DATA1: [%s] [%s] [%s] [%s] [%s]' % (sfKey, tableName, tagName, ccpnData, starInfo)

        emptyFlag = True

        if None in ccpnData:
          ccpnData.remove(None)
        elif '' in ccpnData:
          ccpnData.remove('')
        else:
          emptyFlag = False

        if not ccpnData:
          finalInfoDict[sfKey + '_' + tableName + '_' + tagName] = (STAR_TO_CCPN_TABLE_NAMES[tableName] + ' [' + tagName + ']', 'No "' + starInfo + '" attribute set/selected', orangeColour, False, 'Table Attr', 1, sfKey, tableName, tagName)
        elif emptyFlag:
          finalInfoDict[sfKey + '_' + tableName + '_' + tagName] = (STAR_TO_CCPN_TABLE_NAMES[tableName] + ' [' + tagName + ']', 'Some "' + starInfo + '" attributes not set/selected', yellowColour, False, 'Table Attr', 2, sfKey, tableName, tagName)
        #else:
        #  finalInfoDict[sfKey + '_' + tableName + '_' + tagName] = (sfKey + '/' + tableName + ': ' + tagName, starInfo + ' attribute set/selected', greenColour, True, 'Table Attr')

  finalInfoList = []

  for key in sorted(finalInfoDict.keys() ):
    finalInfoList.append(finalInfoDict[key])

  finalInfoList.sort(listCol6Cmp)

  if not textReport:
    finalInfoList2 = []

    for row in finalInfoList:
      finalInfoList2.append(row[:3])

    return finalInfoList2

  else:
    finalInfoList3 = []

    for row in finalInfoList:
      finalInfoList3.append(row)

    return finalInfoList3

def listCol6Cmp(x, y):
  return cmp(x[5], y[5])
  

if __name__ == '__main__':

  argv = sys.argv[:]
  argc = len(argv)

  mr = None

  if argc > 1:
    projectDir = argv[1]
    mr = loadProject(projectDir)
  else:
    projectDir = None

  nmrEntry = mr.currentNmrEntryStore.findFirstEntry()

  if nmrEntry:

    finalInfoList = checkNmrEntryCompleteness(nmrEntry)

    #print finalInfoList

    for row in finalInfoList:
      print row #[0], row[3:]
