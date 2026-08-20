import os, re, sys

from memops.general.Util import returnMemopsLine, returnMemopsText #, returnMemopsWord
from memops.universal.Util import returnFloat, returnInt

from memops.api import Implementation
from memops.api.Implementation import MemopsRoot

from ccp.api.general import Citation

from ccp.format.pdb import coordinatesIO
from ccpnmr.format.converters.PdbFormat import PdbFormat
from ccpnmr.format.process.matchResonToMolSys import matchCoordAtomsToMolSys

from ccp.general.Io import getChemCompArchiveDataDir
chemCompArchiveDataPath = getChemCompArchiveDataDir()

DEBUG = False

class ReadPdb:

  valueUnitPatt = re.compile("(\d+\.?\d*)\s*([A-Za-z]+)")
  chainCodePatt = re.compile("([Cc][Hh][Aa][Ii][Nn])?\:?\s+([A-Za-z])")

  softVersPatt = re.compile("^(.+?)\s([\d\.]+)$")

  multiSpacePatt = re.compile("\s{2,}")
  emptyLinePatt = re.compile("^(\s*)$")

  lineWithColonPatt = re.compile("^.*:.*")
  lineWithSemiColonPatt = re.compile("^.*;.*")

  concPatt = re.compile("([0-9\.]+)\s*(.*)")

  fieldPatt = re.compile("([0-9\.]+)\s*(.*)")

  unpPatt = re.compile("^[A-Z0-9]{5,6}$") # TODO - what is this in reality
  pdbPatt = re.compile("^[0-9][A-Z0-9]{3}$")

  mutationPatt = re.compile("YES \[(\S+)]")

  remarkTypeToId = {'REFINEMENT.': '3',
                    'EXPERIMENTAL DETAILS': '210'}

  refinePatts = {'program': re.compile("^PROGRAM\s+: (.+)$"),
                 'authors': re.compile("^AUTHORS\s+: (.+)$"),
                 'remarks': re.compile("^OTHER REFINEMENT REMARKS: (.+)$")}

  expPatts = {'expType':   re.compile("^EXPERIMENT TYPE\s+: (.+)$"),
              'expTemp':   re.compile("^TEMPERATURE\s+\(KELVIN\) : (.+)$"),
              'expPh':     re.compile("^PH\s+: (.+)$"),
              'expIon':    re.compile("^IONIC STRENGTH\s+: (.+)$"),
              'expPress':  re.compile("^PRESSURE\s+: (.+)$"),
              'expSample': re.compile("^SAMPLE CONTENTS\s+: (.+)$"),
              'nmrExpts':  re.compile("^NMR EXPERIMENTS CONDUCTED\s+: (.+)$"),
              'specField': re.compile("^SPECTROMETER FIELD STRENGTH\s+: (.+)$"),
              'specModel': re.compile("^SPECTROMETER MODEL\s+: (.+)$"),
              'specManu':  re.compile("^SPECTROMETER MANUFACTURER\s+: (.+)$"),
              'strucInfo': re.compile("^STRUCTURE DETERMINATION.\s*$"),
              'softType':  re.compile("^SOFTWARE USED\s+: (.+)$"),
              'methType':  re.compile("^METHOD USED\s+: (.+)$"),
              'confCalc':  re.compile("^CONFORMERS, NUMBER CALCULATED\s+: (.+)$"),
              'confSub':   re.compile("^CONFORMERS, NUMBER SUBMITTED\s+: (.+)$"),
              'confSel':   re.compile("^CONFORMERS, SELECTION CRITERIA\s+: (.+)$"),
              'confBest':  re.compile("^BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE\s+: (.+)$"),
              'expRemark': re.compile("^REMARK: (.+)$")}

  expPattsList = ['expType', 'expTemp', 'expPh', 'expIon', 'expPress', 'expSample',
                  'nmrExpts', 'specField', 'specModel', 'specManu', 'strucInfo',
                  'softType', 'methType', 'confCalc', 'confSub', 'confSel', 'confBest', 'expRemark']

  expMultiPattsList = ['expSample', 'expTemp', 'expPh', 'expIon', 'expPress',
                       'nmrExpts', 'specField', 'specModel', 'specManu']

  specMultiPattsList = ['specField', 'nmrExpts'] #, 'specModel', 'specManu']


  def __init__(self, pdbFile, ccpnProject, name='readPdb', tkObj=None):

    self.pdbFile = pdbFile
    self.format = 'nmrStar'

    self.ccpnProject = ccpnProject

    self.entryName = name

    self.setCcpnData()

    self.readPdb(tkObj)

    self.setEntryMolSystem()

    self.convert()


  def setCcpnData(self):

    if self.ccpnProject is None:
      self.setCcpnProject()

    self.setNmrProject()

    self.setStructureGeneration()

    self.setCcpnEntry()

    self.setCcpnStores()

    self.setCcpnMolSystem()


  def setCcpnProject(self):

    self.ccpnProject = MemopsRoot(name=self.entryName) # NEW

    if DEBUG:
      print 'NEW PROJECT: [%s]' % self.ccpnProject


  def setNmrProject(self):

    #self.nmrProject = self.ccpnProject.newNmrProject(name=self.entryName) # NEW

    #if DEBUG:
    #  print 'NEW NMR_PROJECT: [%s]' % self.nmrProject

    self.nmrProject = self.ccpnProject.currentNmrProject = self.getTopLevelStore(self.ccpnProject,
                                                                                 'NmrProject',
                                                                                 defaultName=self.entryName)


  def setStructureGeneration(self):

    self.nmrConstraintStore = self.ccpnProject.findFirstNmrConstraintStore(nmrProject=self.nmrProject)

    if not self.nmrConstraintStore:
      self.nmrConstraintStore = self.ccpnProject.newNmrConstraintStore(nmrProject=self.nmrProject) # NEW

      if DEBUG:
        print 'NEW NMR_CONST_STORE: [%s]' % self.nmrConstraintStore

    self.strucGen = self.nmrProject.findFirstStructureGeneration()

    if not self.strucGen:
      sgTitle = 'structure_generation_new_1'

      self.strucGen = self.nmrProject.newStructureGeneration(name=sgTitle,
                                                             generationType='refinement',
                                                             nmrConstraintStore=self.nmrConstraintStore) # NEW

      if DEBUG:
        print 'NEW STRUC_GEN: [%s]' % self.strucGen


  def setCcpnEntry(self):

    self.nmrEntryStore = self.ccpnProject.newNmrEntryStore(name=self.entryName) # NEW

    if DEBUG:
      print 'NEW NMR_ENTRY_STORE: [%s]' % self.nmrEntryStore

    self.ccpnEntry = self.nmrEntryStore.newEntry(name=self.entryName) # NEW

    if DEBUG:
      print 'NEW NMR_ENTRY: [%s]' % self.ccpnEntry


  def setCcpnStores(self):

    self.ccpnProject.currentAffiliationStore = self.getTopLevelStore(self.ccpnProject,
                                                                     'AffiliationStore',
                                                                     defaultName=self.entryName)

    self.ccpnProject.currentCitationStore = self.getTopLevelStore(self.ccpnProject,
                                                                  'CitationStore',
                                                                  defaultName=self.entryName)

    self.ccpnProject.currentTaxonomy = self.getTopLevelStore(self.ccpnProject,
                                                             'Taxonomy',
                                                             defaultName=self.entryName)

    self.ccpnProject.currentSampleStore = self.getTopLevelStore(self.ccpnProject,
                                                                'SampleStore',
                                                                defaultName=self.entryName)

    self.ccpnProject.currentMethodStore = self.getTopLevelStore(self.ccpnProject,
                                                                'MethodStore',
                                                                defaultName=self.entryName)

    self.ccpnProject.currentClassification = self.getTopLevelNSStore(self.ccpnProject,
                                                                     'Classification',
                                                                     defaultNamingSystem=self.entryName)

    self.ccpnProject.currentInstrumentStore = self.getTopLevelStore(self.ccpnProject,
                                                                    'InstrumentStore',
                                                                    defaultName=self.entryName)

    self.ccpnProject.currentRefSampleComponentStore = self.getTopLevelStore(self.ccpnProject,
                                                                            'RefSampleComponentStore',
                                                                            defaultName=self.entryName)

    self.ccpnProject.currentNameMappingStore = self.getTopLevelStore(self.ccpnProject,
                                                                     'NameMappingStore',
                                                                     defaultName=self.entryName)

    self.ccpnProject.currentDatabase = self.getTopLevelStore(self.ccpnProject,
                                                             'Database',
                                                             defaultName=self.entryName)


  def setCcpnMolSystem(self):

    self.molSystem = self.ccpnProject.newMolSystem(code=self.entryName) # NEW

    if DEBUG:
      print 'NEW MOL_SYSTEM: [%s]' % self.molSystem


  def readPdb(self, tkObj=None):

    self.pdbDataFile = None

    self.pdbFormat = PdbFormat(self.ccpnProject, tkObj, allowPopups=False)

    from ccp.format.pdb.coordinatesIO import PdbCoordinateFile

    pdbFile = PdbCoordinateFile(self.pdbFile)
    pdbFile.read(maxNum=1)

    # TODO: How do we make it deal with complex heteromeric complexes?

    numIdenticalChains = len(pdbFile.molIdToChains[1])

    self.pdbTitle = pdbFile.title.strip()

    if self.multiSpacePatt.search(self.pdbTitle):
      (self.pdbTitle, _numSubs) = self.multiSpacePatt.subn(' ', self.pdbTitle)

    self.ccpnEntry.title = self.pdbTitle
    self.molSystem.name = self.pdbTitle

    self.pdbFormat.readSequence(self.pdbFile,
                                molSystem=self.molSystem,
                                verbose=True,
                                minimalPrompts=True,
                                chemCompPath=chemCompArchiveDataPath,
                                numIdenticalChains=numIdenticalChains) # NEW

    if DEBUG:
      print 'INP MOL_SYSTEM: [%s]' % self.molSystem

    self.pdbFormat.readPeopleAndCitations(self.pdbFile, doNotMakeCcpnObjects=True, verbose=True, minimalPrompts=True)

    if DEBUG:
      print 'INP PEOP_CIT (PDB_FORMAT): [%s]' % self.pdbFormat

    self.pdbDataFile = self.pdbFormat.peopleAndCitationsFile.pdbFile


  def setEntryMolSystem(self):

    self.ccpnEntry.molSystem = self.molSystem # SET

    if DEBUG:
      print 'SET MOL_SYSTEM: [%s] [%s]' % (self.ccpnEntry.molSystem, self.ccpnEntry)


  def convert(self):

    self.setTitleAndKeywords()

    self.setMolInfo()

    self.setMolecules()

    self.setDbRefs()

    self.setStructInfo()

    self.setPeopleAndCitationInfo()

    self.setNmrExperiments()

    self.setMolComponents()


  def setTitleAndKeywords(self):

    self.ccpnEntry.setTitle(self.pdbTitle) # CHANGE

    if DEBUG:
      print 'SET TITLE: [%s] [%s]' % (self.ccpnEntry, self.pdbTitle)

    if self.pdbDataFile.headerVars.has_key('Keywds'):
      kywds = [ kwd.strip() for kwd in self.pdbDataFile.headerVars['Keywds'].split(',') ]

      for kywd in kywds:
        self.ccpnEntry.addKeyword(kywd) # ADD

        if DEBUG:
          print 'ADD KEYWD: [%s] [%s]' % (self.ccpnEntry, kywd)


  def setMolInfo(self):

    chains = self.pdbDataFile.chains

    sourceInfo = self.pdbDataFile.sourceInfo

    defaultChainCodes = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    seenChains = []

    for chainNum in range(len(chains) ):

      chainId = chains[chainNum-1].molId
      chainCodeFind = ''

      chainCodeText = chains[chainNum-1].chainId

      searchObj = self.chainCodePatt.search(chainCodeText)

      if searchObj:
        chainCode = searchObj.group(2)

      elif chainCodeText in defaultChainCodes:
        chainCode = chainCodeText

      else:
        chainCode = defaultChainCodes[(int(chainId) - 1)]

      chain = self.molSystem.findFirstChain(code=chainCode)

      if not chain:
        print "Error: no chain found for code '%s'. No source/expression information set." % chainCode
        continue

      chainCodeFind = chainCode

      if chain not in seenChains:
        seenChains.append(chain)
      else:
        continue

      molecule = chain.molecule

      pdbChain = chains[chainNum-1]

      if hasattr(pdbChain, 'ecNums'):
        ecNumStr = ','.join(pdbChain.ecNums)
        if ecNumStr != '':
          msn = molecule.newMoleculeSysName(name=ecNumStr, namingSystem='EC') # NEW

          if DEBUG:
            print 'NEW: MOL_SYS_NAME: [%s] [%s]' % (molecule, msn)

      if hasattr(pdbChain, 'mutation') and pdbChain.mutation:
        searchObj = self.mutationPatt.search(pdbChain.mutation)

        mutation = pdbChain.mutation

        if searchObj:
          mutation = searchObj.group(1)

        keywds = {'application': self.format,
                  'keyword':     'Mutation',
                  'value':       mutation}

        appData = Implementation.AppDataString(**keywds)
        molecule.addApplicationData(appData)

        if DEBUG:
          print 'ADD MUT: [%s] [%s]' % (molecule, mutation)

      if hasattr(pdbChain, 'fragment') and pdbChain.fragment:
        keywds = {'application': self.format,
                  'keyword':     'Fragment',
                  'value':       pdbChain.fragment}

        appData = Implementation.AppDataString(**keywds)
        molecule.addApplicationData(appData)

        if DEBUG:
          print 'ADD FRAG: [%s] [%s]' % (molecule, pdbChain.fragment)

      if len(sourceInfo) >= chainNum:

        keywds = {
          'organismName':   sourceInfo[chainNum-1].get('ORGANISM_COMMON', 'n/a'),
          'scientificName': sourceInfo[chainNum-1].get('ORGANISM_SCIENTIFIC', 'n/a'),
          'strain':         sourceInfo[chainNum-1].get('STRAIN'),
          'cellLine':       sourceInfo[chainNum-1].get('CELL_LINE'),
          'organ':          sourceInfo[chainNum-1].get('ORGAN'),
          'variant':        sourceInfo[chainNum-1].get('VARIANT'),
          'tissue':         sourceInfo[chainNum-1].get('TISSUE'),
          'details':        sourceInfo[chainNum-1].get('OTHER_DETAILS'),
          'atccNumber':     sourceInfo[chainNum-1].get('ATCC'),
          'cellType':       sourceInfo[chainNum-1].get('CELL'),
          'organelle':      sourceInfo[chainNum-1].get('ORGANELLE'),
          'geneMnemonic':   sourceInfo[chainNum-1].get('GENE'),
          'ncbiTaxonomyId': sourceInfo[chainNum-1].get('ORGANISM_TAXID'),
          }

      else:

        keywds = {'organismName':   'n/a',
                  'scientificName': 'n/a'}

      naturalSource = self.ccpnProject.currentTaxonomy.newNaturalSource(**keywds) # NEW

      if DEBUG:
        print 'NEW NAT_SRC: [%s]' % naturalSource

      molecule.setNaturalSource(naturalSource) # SET

      if DEBUG:
        print 'SET NAT_SRC: [%s] [%s]' % (molecule, naturalSource)


      syn = sourceInfo[chainNum-1].get('SYNTHETIC')

      if syn:
        prodMeth = 'chemical synthesis'
      else:
        prodMeth = 'recombinant technology'

      entryMol = self.ccpnEntry.findFirstEntryMolecule(molecule=molecule, productionMethod=prodMeth)

      if not entryMol:
        entryMol = self.ccpnEntry.newEntryMolecule(molecule=molecule, productionMethod=prodMeth) # NEW

        if DEBUG:
          print 'NEW ENTRY_MOL: [%s] [%s]' % (self.ccpnEntry, entryMol)

      experimentalSource = None

      if len(sourceInfo) >= 1 and sourceInfo[chainNum-1]:
        vectorType = sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_VECTOR')

        if vectorType and vectorType != '':
          entryMol.vectorType = vectorType.strip().lower()

        keywds2 = {
          #'organismName':   sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_COMMON', 'unknown'),
          'scientificName': sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM', 'unknown'),
          'strain':         sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_STRAIN'),
          'cellLine':       sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_CELL_LINE'),
          'organ':          sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_ORGAN'),
          'variant':        sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_VARIANT'),
          'tissue':         sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_TISSUE'),
          'details':        sourceInfo[chainNum-1].get('OTHER_DETAILS'),
          'plasmid':        sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_PLASMID'),
          'atccNumber':     sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_ATCC'),
          'cellType':       sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_CELL'),
          'geneMnemonic':   sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_GENE'),
          'ncbiTaxonomyId': sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_TAXID'),
          }

        comName = sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM_COMMON')

        if not comName:
          comName = sourceInfo[chainNum-1].get('EXPRESSION_SYSTEM', 'unknown')

        keywds2['organismName'] = comName


      else:

        keywds2 = {'organismName':   'unknown',
                   'scientificName': 'unknown'}

      experimentalSource = self.ccpnProject.currentTaxonomy.findFirstNaturalSource(**keywds2)

      if not experimentalSource:
        experimentalSource = self.ccpnProject.currentTaxonomy.newNaturalSource(**keywds2) # NEW

        if DEBUG:
          print 'NEW EXP_SRC: [%s]' % experimentalSource

      entryMol.setExperimentalSource(experimentalSource) # SET

      if DEBUG:
        print 'SET EXP_SRC: [%s] [%s]' % (entryMol, experimentalSource)


  def setMolecules(self):

    for chain in self.molSystem.sortedChains():
      molecule = chain.molecule
      molecule.isFinalised = True


  def setDbRefs(self):

    for dbRef in self.pdbDataFile.dbRefs:
      #print 'DBREF: [%s]' % str(dbRef)

      dbRefList = [ dbr.strip() for dbr in dbRef ]

      dbCode = dbRefList[7]

      dbId = ''
      if dbRefList[8][0:3] != 'XXX':
        dbId = dbRefList[8]

      database = None

      searchObj = self.unpPatt.search(dbCode)

      if searchObj:
        database = self.ccpnProject.findFirstDatabase(name='UniProt')

        if not database:
          database = self.ccpnProject.findFirstDatabase(name='UNP')

      if not database:
        searchObj = self.pdbPatt.search(dbCode)

        if searchObj:
          database = self.ccpnProject.findFirstDatabase(name='PDB')

      origDbName = dbName = dbRefList[6]

      if dbName in ('UNP', 'SWS', 'UniProt'):
        dbName = 'UniProt'

      if database:
        if database.name != dbName:
          database = self.ccpnProject.findFirstDatabase(name=dbName)
          if not database:
            database = self.ccpnProject.newDatabase(name=dbName) # NEW

            if DEBUG:
              print 'NEW DATABASE: [%s]' % database

      if not database and dbName and dbName != '' and dbName.lower() != 'n/a' and dbName[0:3] != 'XXX':
        database = self.ccpnProject.findFirstDatabase(name=dbName)
        if not database:
          database = self.ccpnProject.newDatabase(name=dbName)

          if DEBUG:
            print 'NEW DATABASE: [%s]' % database

      if not database:
        database = self.ccpnProject.currentDatabase

        if DEBUG:
          print 'CUR DATABASE: [%s]' % database

      if dbCode != 'n/a' and database:
        keywds = {}

        myEntry = database.newEntry(name=dbCode) # NEW

        if DEBUG:
          print 'NEW DB_ENTRY: [%s] [%s]' % (database, myEntry)

          if dbId != '':
            myEntry.code = dbId
          else:
            myEntry.code = dbCode

        startPdb = dbRefList[2]
        if startPdb and not startPdb.startswith('X'):
          keywds['dbRefAlignBegin'] = returnInt(startPdb)

        endPdb = dbRefList[4]
        if endPdb and not endPdb.startswith('X'):
          keywds['dbRefAlignEnd'] = returnInt(endPdb)

        startRef = dbRefList[9]
        endRef   = dbRefList[11]

        if startRef or endRef:
          keywds['details'] = startRef + ':' + endRef

        chainId = dbRefList[1]

        chain = self.molSystem.findFirstChain(code=chainId)

        if not chain:
          chain = self.molSystem.findFirstChain(code=' ')

        if chain:
          molecule = chain.molecule

          alignment = molecule.newAlignment(dbRef=myEntry, **keywds) # NEW

          if DEBUG:
            print 'NEW ALN: [%s]' % alignment


  def setStructInfo(self):

    coordinateFile = coordinatesIO.PdbCoordinateFile(self.pdbFile)
    coordinateFile.read(maxNum=1)

    modelCoordKeys = coordinateFile.modelCoordinates.keys()
    modelCoordKeys.sort()

    refModel = modelCoordKeys[0]
    coords = coordinateFile.modelCoordinates[refModel]

    keywds = {}

    print "Trying to automap coordinate atoms."

    forceChainMappings = matchCoordAtomsToMolSys(coords, self.molSystem, test=0)

    if forceChainMappings:
      keywds['forceChainMappings'] = forceChainMappings

      keywds['ignoreUnknownChemComps'] = True
      keywds['forceReadSequence'] = False
      keywds['linkAtoms'] = False
      keywds['minimalPrompts'] = 1
      keywds['chemCompPath'] = chemCompArchiveDataPath
      keywds['maxNum'] = 999

    self.pdbFormat.readCoordinates(self.pdbFile, molSystem=self.molSystem, strucGen=self.strucGen, **keywds) # NEW

    if DEBUG:
      print 'INP STRUC_GEN: [%s]' % self.strucGen

    if not self.ccpnEntry in self.strucGen.sortedEntries():
      self.strucGen.addEntry(self.ccpnEntry) # ADD

      if DEBUG:
        print 'ADD STRUC_GEN: [%s] [%s]' % (self.ccpnEntry, self.strucGen)


  def setPeopleAndCitationInfo(self):

    authors = self.pdbFormat.peopleAndCitationsFile.authors

    for person in authors:

      authorName = ''

      if person.familyName is not None:
        authorName += person.familyName

      if person.firstName is not None:
        authorName = person.firstName + '.' + authorName

      elif person.initials is not None:
        authorName = '.'.join(person.initials) + '.' + authorName

      author = self.getVaguePerson(authorName)

      if not author in self.ccpnEntry.authors:
        self.ccpnEntry.addAuthor(author) # ADD

        if DEBUG:
          print 'ADD AUTHOR: [%s] [%s]' % (self.ccpnEntry, author)

    pdbCitations = self.pdbFormat.peopleAndCitationsFile.citations

    for citation in pdbCitations:

      title = ''

      if citation.title:
        title = citation.title.replace(os.linesep, ' ').strip()

        if self.multiSpacePatt.search(title):
          (title,_numSubs) = self.multiSpacePatt.subn(' ', title)
          title = title.strip()

      keywds = {'title': title}

      keywds['authors'] = []

      authors = citation.authors

      for person in authors:

        authorName = ''

        if person.familyName is not None:
          authorName += person.familyName

        if person.firstName is not None:
          authorName = person.firstName + '.' + authorName

        elif person.initials is not None:
          authorName = '.'.join(person.initials) + '.' + authorName

        author = self.getVaguePerson(authorName)

        if not author in keywds['authors']:
          keywds['authors'].append(author)

      keywds['editors'] = []

      editors = citation.editors

      for person in editors:

        editorName = ''

        if person.familyName is not None:
          editorName += person.familyName

        if person.firstName is not None:
          editorName = person.firstName + '.' + editorName

        elif person.initials is not None:
          editorName = '.'.join(person.initials) + '.' + editorName

        editor = self.getVaguePerson(editorName)

        if not editor in keywds['editors']:
          keywds['editors'].append(editor)

      citType = 'JournalCitation'
      status = citation.status

      if status == 'to be published':
        status = 'in preparation'

      keywds['status'] = status
      keywds['firstPage'] = citation.firstPage
      keywds['lastPage'] = citation.lastPage
      keywds['year'] = citation.year

      if citation.details:
        keywds['details'] = citation.details

      if citation.type == 'journal':

        citType = 'JournalCitation'

        keywds['journalAbbreviation'] = citation.pubShortName
        keywds['journalFullName'] = citation.pubLongName
        if citation.volume != '':
          keywds['volume'] = citation.volume

        if DEBUG:
          print 'SETTING AS journal'

      elif citation.type == 'book':

        citType = 'BookCitation'

        keywds['publisher'] = citation.publisher
        keywds['publisherCity'] = citation.publisherPlace
        if citation.volume != '':
          keywds['volume'] = citation.volume

        if DEBUG:
          print 'SETTING AS book'

      elif citation.type == 'thesis':

        citType = 'ThesisCitation'

        keywds['city'] = citation.publisherPlace
        keywds['institution'] = citation.publisher

        if DEBUG:
          print 'SETTING AS thesis'

      else:

        if DEBUG:
          print 'DEFAULTING to journal'

      ccpnCitation = getattr(Citation,citType)(self.ccpnProject.currentCitationStore, **keywds) # NEW

      if DEBUG:
        print 'NEW CIT: [%s]' % ccpnCitation

      if citation.isPrimary:
        self.ccpnEntry.setPrimaryCitation(ccpnCitation) # SET

        if DEBUG:
          print 'SET PRIM_CIT: [%s] [%s]' % (self.ccpnEntry, ccpnCitation)

      else:
        self.ccpnEntry.addOtherCitation(ccpnCitation) # ADD

        if DEBUG:
          print 'ADD OTHER_CIT: [%s] [%s]' % (self.ccpnEntry, ccpnCitation)


  def setNmrExperiments(self):

    remarks = self.pdbDataFile.remarks

    expInfo = self.parseExpLines(remarks)

    sampleCategoryName = self.entryName

    classification = self.ccpnProject.currentClassification

    sampleCategory = classification.newSampleCategory(name=sampleCategoryName) # NEW

    if DEBUG:
      print 'NEW SAMP_CAT: [%s] [%s]' % (classification, sampleCategory)

    nmrType = 'solutionNMR' # TODO: We can get this from the PDB file I think.
    
    if nmrType == 'solutionNMR':
      sampleState = 'isotropic'
    else:
      sampleState = 'solid'

    multiCondsFlag = False
    multiSpecsFlag = False

    numSamples = 1

    patt = self.expMultiPattsList[0] # Samples
    patt2 = self.expMultiPattsList[1] # Temps

    if self.lineWithSemiColonPatt.search(expInfo[patt]):
      multiCondsFlag = True

      numSamples = expInfo[patt].count(';') + 1
      numSamples2 = expInfo[patt2].count(';') + 1

      #print 'NUM: [%s] [%s]' % (numSamples, numSamples2)

    if numSamples == 1:
      patt = self.specMultiPattsList[0]

      if self.lineWithSemiColonPatt.search(expInfo[patt]):
        multiSpecsFlag = True

        numSamples = expInfo[patt].count(';') + 1

    if multiCondsFlag or multiSpecsFlag:
      for patt in self.expMultiPattsList:
        if expInfo[patt].count(';'):
          expInfo[patt] = [ field.strip() for field in expInfo[patt].split(';') ]
        else:
          expInfo[patt] = [ expInfo[patt] ]

    else:
      for patt in self.expMultiPattsList:
        expInfo[patt] = [ expInfo[patt] ]

    scSeenFlag   = False
    scNewSamFlag = False

    for i in range(numSamples):

      # TODO This 'if' statement is because the PDB file format is wrong at the moment - eventually, we should just need the 'if' branch.

      if len(expInfo['expSample']) > i:
        sampleName = 'sample_new_' + str(i + 1) #expInfo['expSample'][i]
        solventSys = expInfo['expSample'][i]
      else:
        sampleName = 'sample_new_1' #expInfo['expSample'][0]
        solventSys = expInfo['expSample'][0]

      if len(expInfo['expPh']) > i:
        ph = expInfo['expPh'][i]
      else:
        ph = expInfo['expPh'][0]

      if ph and ph != 'NULL':
        ph = returnFloat(ph)
      else:
        ph = 0.0

      if len(expInfo['expIon']) > i:
        searchObj = self.concPatt.search(expInfo['expIon'][i])
      else:
        searchObj = self.concPatt.search(expInfo['expIon'][0])

      ionicStrength = 0.0

      if searchObj:
        ionicStrength = searchObj.group(1).strip()

        if ionicStrength and ionicStrength != 'NULL':
          ionicStrength = returnFloat(ionicStrength)
        else:
          ionicStrength = 0.0

      # Definitely can't parse sample components from PDB files at the moment.

      nmrCompoundConc = solventSys

      keywds = {
        'ph':            ph,
        'ionicStrength': ionicStrength,
        'details':       nmrCompoundConc
        }

      sample = self.ccpnProject.currentSampleStore.findFirstAbstractSample(**keywds) 

      keywds['name'] = sampleName

      if not sample:
        sample = self.ccpnProject.currentSampleStore.newSample(sampleCategories=[sampleCategory], **keywds) # NEW - already in the project?

        if DEBUG:
          print 'NEW SAMPLE: [%s]' % sample

      appData = sample.findFirstApplicationData(application=self.format,keyword='solventSys')

      if not appData and solventSys and solventSys != 'NULL':

        keywds = {'application': self.format,
                  'keyword':     'solventSys',
                  'value':       solventSys}

        appData = Implementation.AppDataString(**keywds)
        sample.addApplicationData(appData)

        if DEBUG:
          print 'ADD SOLV: [%s] [%s]' % (sample, solventSys)

      if len(expInfo['expTemp']) > i:
        temp = expInfo['expTemp'][i]
      else:
        temp = expInfo['expTemp'][0]
        
      if temp and temp != 'NULL':
        pass #temp = returnFloat(temp)
      else:
        temp = None

      if len(expInfo['expPress']) > i:
        pressure = expInfo['expPress'][i]
      else:
        pressure = expInfo['expPress'][0]

      if pressure and pressure != 'NULL':
        #pressure = returnFloat(pressure)
        if pressure == 'AMBIENT':
          pressure = '1.0 atm'

      else:
        pressure = None

      details = "pH [%s], temp [%s], pressure [%s], ionStrength [%s]" % (ph,
                                                                         temp,
                                                                         pressure,
                                                                         ionicStrength)

      sampleConditionSet = self.nmrProject.findFirstSampleConditionSet(details=details)

      if not sampleConditionSet:
        scsName = 'sample_conditions_new' + str(i + 1)

        sampleConditionSet = self.nmrProject.newSampleConditionSet(name=scsName, details=details) # NEW

        if DEBUG:
          print 'NEW SCS: [%s]' % sampleConditionSet

        for (valueName, value, defaultUnit) in ( ('pH', str(ph), 'pH'),
                                                 ('pressure', str(pressure), 'atm'),
                                                 ('temperature', str(temp), 'K'),
                                                 ('ionic strength', str(ionicStrength), 'M') ):

          if not value or value in ('None', 'NULL'):
            continue

          searchObj = self.valueUnitPatt.search(value)

          unit = None

          if searchObj:
            value = searchObj.group(1)

            unit = searchObj.group(2)

          if not unit:
            unit = defaultUnit

          if value:
            if value.count('.'):
              value = returnFloat(value)
            else:
              value = returnInt(value)

          newSamCond = sampleConditionSet.findFirstSampleCondition(condition=valueName, value=value, unit=unit)

          if not newSamCond:
            newSamCond = sampleConditionSet.newSampleCondition(condition=valueName, value=value, unit=unit) # NEW

            if DEBUG:
              print 'NEW SC: [%s] [%s]' % (sampleConditionSet, newSamCond)

        """
        for scs in self.nmrProject.sortedSampleConditionSets():
          if scs != sampleConditionSet:

            for sc in scs.sortedSampleConditions():
              if sampleConditionSet.findFirstSampleCondition(condition=sc.condition, value=sc.value, unit=sc.unit):
                scSeenFlag = True

              elif sc.condition == 'ionic strength':
                newSc = sampleConditionSet.findFirstSampleCondition(condition=sc.condition)

                if not newSc:
                  scSeenFlag = True # ???
                  break

                if ( (sc.unit == 'M' and newSc.unit == 'mM' and (sc.value*1000) == newSc.value) or 
                        (sc.unit == 'mM' and newSc.unit == 'M' and (newSc.value*1000) == sc.value) ):
                  scSeenFlag = True

              elif sc.condition not in ('pH', 'pressure', 'temperature', 'ionic strength'):
                scSeenFlag = True # ???
                break

              else:
                scSeenFlag = False
                break

            if scSeenFlag:
              sampleConditionSet = scs
              break
        """

      if len(expInfo['specManu']) > i:
        manufacturerName = expInfo['specManu'][i]
        specModel = expInfo['specModel'][i]
        specField = expInfo['specField'][i]

      else:
        pressure = expInfo['specManu'][0]
        specModel = expInfo['specModel'][0]
        specField = expInfo['specField'][0]

      manufacturer = self.ccpnProject.currentAffiliationStore.findFirstOrganisation(name=manufacturerName)

      if not manufacturer:
        manufacturer = self.ccpnProject.currentAffiliationStore.newOrganisation(name=manufacturerName) # NEW

        if DEBUG:
          print 'NEW SPEC_MANU: [%s]' % manufacturer

      searchObj = self.fieldPatt.search(str(specField) )

      if searchObj:
        specField = searchObj.group(1)

        if not specField or specField == 'NULL':
          specField = None

      specName = "%s %s-%s" % (manufacturerName, specModel, specField)

      keywds = {
        'nominalFreq':  specField,
        'protonFreq':   returnFloat(specField), # Can't get exact field
        'name':         specName,
        'model':        specModel,
        'manufacturer': manufacturer
        }

      if not specName or specName == '':
        specName = 'spectrometer_new_' + str(i + 1)

      keywds['name'] = specName

      spectrometer = self.ccpnProject.currentInstrumentStore.findFirstInstrument(className='NmrSpectrometer', **keywds)

      if not spectrometer:
        spectrometer = self.ccpnProject.currentInstrumentStore.newNmrSpectrometer(**keywds) # NEW

        if DEBUG:
          print 'NEW SPEC: [%s]' % spectrometer

      if len(expInfo['nmrExpts']) <= i:
        continue

      if len(expInfo['nmrExpts']) > i:
        nmrExpString = expInfo['nmrExpts'][i]
      else:
        nmrExpString = expInfo['nmrExpts'][0]

      if nmrExpString:
        nmrExpString = nmrExpString.replace(os.linesep, '')

        if nmrExpString != 'NULL':
          for separator in (',',';','and','AND'):
          #for separator in (';','and','AND'):
            nmrExpStrings = nmrExpString.split(separator)

            if len(nmrExpStrings) > 1:
              break

          if len(nmrExpStrings) == 1:
            nmrExpStrings = [nmrExpString]

        else:
          nmrExpStrings = []

      exptNames = []

      for expt in self.nmrProject.sortedExperiments():
        if expt.name.upper() not in exptNames:
          exptNames.append(expt.name.upper() )

      for nmrExpName in nmrExpStrings:

        nmrExpName = nmrExpName.strip()

        if not nmrExpName:
          nmrExpName = 'nmrExpt'
        elif len(nmrExpName) >= 80:
          nmrExpName = nmrExpName[:80]

        scsName = ''
        if sampleConditionSet:
          scsName = sampleConditionSet.name

        nmrSpecName = ''
        if spectrometer:
          nmrSpecName = spectrometer.name

        keywds = {
          'numDim':       2, # Cannot get this info from names.
          'name':         nmrExpName,
          'spectrometer': spectrometer,
          'sampleState':  sampleState # Defined higher up.
          }

        nmrExp = None

        if nmrExpName.upper() in exptNames:
          for expt in self.nmrProject.sortedExperiments():
            if nmrExpName.upper() == expt.name.upper():
              expSamName = ''
              if expt.sample:
                expSamName = expt.sample.name

                expScsName = ''
              if expt.sampleConditionSet:
                expScsName = expt.sampleConditionSet.name

                expSpecName = ''
              if expt.spectrometer:
                expSpecName = expt.spectrometer.name

              if sampleName == expSamName and scsName == expScsName and specName == expSpecName:
                nmrExp = expt
                break

        if not nmrExp:
          nmrExp = self.nmrProject.newExperiment(**keywds) # NEW

          if DEBUG:
            print 'NEW EXPT: [%s]' % nmrExp

          scNewSamFlag = True

        if not self.molSystem in nmrExp.sortedMolSystems():
          nmrExp.addMolSystem(self.molSystem) # ADD

          if DEBUG:
            print 'ADD MOL_SYS: [%s] [%s]' % (nmrExp, self.molSystem)

        if not self.ccpnEntry in nmrExp.sortedEntries():
          nmrExp.addEntry(self.ccpnEntry) # ADD

          if DEBUG:
            print 'ADD NMR_ENTRY: [%s] [%s]' % (nmrExp, self.ccpnEntry)

        if sample and scNewSamFlag and not nmrExp in sample.sortedNmrExperiments():
          sample.addNmrExperiment(nmrExp) # ADD

          if DEBUG:
            print 'ADD SAMPLE: [%s] [%s]' % (nmrExp, sample)

        if sampleConditionSet and not nmrExp in sampleConditionSet.sortedExperiments():
          sampleConditionSet.addExperiment(nmrExp) # ADD

          if DEBUG:
            print 'ADD SCS: [%s] [%s]' % (nmrExp, sampleConditionSet)

        scsName = ''
        if nmrExp.sampleConditionSet:
          scsName = nmrExp.sampleConditionSet.name

        nmrSpecName = ''
        if nmrExp.spectrometer:
          nmrSpecNameName = nmrExp.spectrometer.name

    refineInfo = self.parseRefineLines(remarks)

    if self.strucGen:
      if not self.strucGen.details and 'expRemark' in expInfo:
        self.strucGen.details = expInfo['expRemark']

      strEns = self.strucGen.structureEnsemble

      if not self.ccpnEntry in self.strucGen.sortedEntries():
        self.strucGen.addEntry(self.ccpnEntry) # ADD

        if DEBUG:
          print 'ADD STRUC_GEN: [%s] [%s]' % (self.ccpnEntry, self.strucGen)

      confBest = None

      if 'confBest' in expInfo:
        confBest = expInfo['confBest']

      if confBest and confBest != 'NULL':

        keywds = {'application': self.format,
                  'keyword':     'representative',
                  'value':       confBest}

        appData = Implementation.AppDataString(**keywds)
        strEns.addApplicationData(appData)

        if DEBUG:
          print 'ADD REPR: [%s] [%s]' % (strEns, confBest)

      confSel = None

      if 'confSel' in expInfo:
        confSel = expInfo['confSel']

      if confSel and confSel != 'NULL':

        keywds = {'application': self.format,
                  'keyword':     'criteria',
                  'value':       confSel}

        appData = Implementation.AppDataString(**keywds)
        strEns.addApplicationData(appData)

        if DEBUG:
          print 'ADD CRITERIA: [%s] [%s]' % (strEns, confSel)

      confCalc = None

      if 'confCalc' in expInfo:
        confCalc = expInfo['confCalc']

      if confCalc and confCalc != 'NULL':

        keywds = {'application': self.format,
                  'keyword':     'calculated',
                  'value':       confCalc}

        appData = Implementation.AppDataString(**keywds)
        strEns.addApplicationData(appData)

        if DEBUG:
          print 'ADD CALC: [%s] [%s]' % (strEns, confCalc)

    softwareList = expInfo['softType']
    methodName = expInfo['methType']
    methodName = methodName and methodName[:79]

    for separator in (',',';','and','AND'):
      softwareStrings = softwareList.split(separator)

      if len(softwareStrings) > 1:
        break

    for i in range(len(softwareStrings) ):
      softwareStrings[i] = softwareStrings[i].strip()

    if 'program' in refineInfo:
      strucCalcSoft = refineInfo['program']

      if strucCalcSoft and strucCalcSoft != 'NULL':
        if strucCalcSoft in softwareStrings:
          softwareStrings.remove(strucCalcSoft)

        authors = None

        if 'authors' in refineInfo and refineInfo['authors'] != 'NULL':
          authors = refineInfo['authors']

        remarks = None

        if 'remarks' in refineInfo and refineInfo['remarks'] != 'NULL':
          remarks = refineInfo['remarks']

        (software, method) = self.getSoftwareMethod(strucCalcSoft, methodName, authors, remarks)

        if software not in self.ccpnEntry.sortedSoftware():
          self.ccpnEntry.addSoftware(software) # ADD

          if DEBUG:
            print 'ADD SOFTWARE: [%s] [%s]' % (self.ccpnEntry, software)

        if self.strucGen and not self.strucGen.method and method:
          self.strucGen.setMethod(method) # SET

          if DEBUG:
            print 'SET METHOD: [%s] [%s]' % (self.strucGen, method)

    for i in range(len(softwareStrings) ):
      softwareName = softwareStrings[i]

      if softwareName and softwareName != 'NULL':
        (software,method) = self.getSoftwareMethod(softwareName)

        if software not in self.ccpnEntry.sortedSoftware():
          self.ccpnEntry.addSoftware(software) # ADD

          if DEBUG:
            print 'ADD SOFTWARE: [%s] [%s]' % (self.ccpnEntry, software)


  def setMolComponents(self):

    seenMols = []

    for chain in self.molSystem.sortedChains():
      
      molecule = chain.molecule

      if molecule in seenMols:
        continue

      else:
        seenMols.append(molecule)

      if len(molecule.molResidues) > 1 and molecule.isStdLinear:
        molName = molecule.name

        for sample in self.ccpnProject.currentSampleStore.sortedAbstractSamples():

          molComp = self.ccpnProject.currentRefSampleComponentStore.findFirstComponent(name=molName + ':' + sample.name)

          if not molComp:
            molComp = self.ccpnProject.currentRefSampleComponentStore.findFirstComponent(name=molName)

          if not molComp:
            molComp = self.ccpnProject.currentRefSampleComponentStore.findFirstComponent(details=molName)

          if not molComp:
            molComp = self.ccpnProject.currentRefSampleComponentStore.newMolComponent(name=molName + ':' + sample.name, details=molName) # NEW

            if DEBUG:
              print 'NEW MOL_COMP: [%s]' % molComp

          natSrc = None

          if molecule.naturalSource:
            natSrc = molecule.naturalSource

          if not molComp.naturalSource and natSrc and molComp.naturalSource != natSrc:
            molComp.setNaturalSource(natSrc) # SET

            if DEBUG:
              print 'SET NAT_SRC: [%s] [%s]' % (molComp, natSrc)

          if molComp not in molecule.sortedMolComponents():
            molecule.addMolComponent(molComp) # ADD

            if DEBUG:
              print 'ADD MOL MOL_COMP: [%s] [%s]' % (molecule, molComp)


  def parseRefineLines(self,remarks):

    refineInfo = {}

    refineString = 'REFINEMENT.'

    if remarks.has_key(self.remarkTypeToId[refineString]):
      refine = remarks[self.remarkTypeToId[refineString] ]

      #print 'REFINE: [%s]' % refine

      refineLines = refine.split(os.linesep)

      lineNum = 0

      if refineLines[lineNum] != refineString:
        print "  Error: No refinement information in header"

      lineNum += 1

      searchObj = self.refinePatts['program'].search(refineLines[lineNum])

      if not refineInfo.has_key('program'):
        refineInfo['program'] = ''

      if searchObj:
        refineInfo['program'] = searchObj.group(1)

      while refineLines[lineNum]:
        lineNum += 1

        line = refineLines[lineNum]

        searchObj = self.refinePatts['authors'].search(line)

        if searchObj:
          break

        refineInfo['program'] += ' ' + line.strip()

      if not refineInfo.has_key('authors'):
        refineInfo['authors'] = ''

      if searchObj:
        refineInfo['authors'] = searchObj.group(1)

      lineNum += 1

      searchObj = None

      for line in refineLines[lineNum:]:
        lineNum += 1

        searchObj = self.refinePatts['remarks'].search(line)
        searchObj2 = self.lineWithColonPatt.search(line)

        line = line.strip()

        #print 'LINE: [%s] [%s] [%s]' % (line, searchObj, searchObj2)

        if searchObj or searchObj2 or line == '' or line[0:9] == 'NUMBER OF':
          break
        else:
          refineInfo['authors'] += line

      if not refineInfo.has_key('remarks'):
        refineInfo['remarks'] = ''

      if not searchObj:
        while refineLines[lineNum]:
          lineNum += 1

          line = refineLines[lineNum]

          searchObj = self.refinePatts['remarks'].search(line)

          if searchObj:
            break

      if searchObj:
        refineInfo['remarks'] = searchObj.group(1)

      lineNum += 1

      for line in refineLines[lineNum:]:
        lineNum += 1

        searchObj = self.lineWithColonPatt.search(line)

        line = line.strip()

        #print 'LINE: [%s] [%s] [%s]' % (line, searchObj, searchObj2)

        if searchObj or line == '':
          break
        else:
          refineInfo['remarks'] += ' ' + line

      #for k in refineInfo.keys():
      #  print 'DATA1: [' + k + '] [' + refineInfo[k] + ']'
      #print

    return refineInfo


  def parseExpLines(self,remarks):

    expInfo = {}

    expString = 'EXPERIMENTAL DETAILS'

    if remarks.has_key(self.remarkTypeToId[expString]):
      expDetails = remarks[self.remarkTypeToId[expString] ]

      expDetLines = expDetails.split(os.linesep)

      #print 'EXP DETAILS: [' + str(expDetLines) + ']'

      lineNum = 0
      fieldCount = 0

      if expDetLines[lineNum] != expString:
        print "  Error: No experimental details in header"

      lineNum += 1

      patt = self.expPattsList[fieldCount]
      line = expDetLines[lineNum]

      searchObj = self.expPatts[patt].search(line)

      if not expInfo.has_key(patt):
        expInfo[patt] = ''

      if searchObj:
        expInfo[patt] = searchObj.group(1)

      lineNum += 1

      for fieldCount in range(1, (len(self.expPattsList) ) ):

        patt = self.expPattsList[fieldCount]

        if not expInfo.has_key(patt):
          expInfo[patt] = ''

        while expDetLines[lineNum]:
          line = expDetLines[lineNum]

          searchObj = self.expPatts[patt].search(line)

          lineNum += 1

          if searchObj:
            break

          lastPatt = self.expPattsList[fieldCount-1]

          line = line.strip()

          if expInfo[lastPatt][-1] in ('-', '.', '(', ')', '/', '\\', '%'):
            expInfo[lastPatt] += '' + line
          else:
            expInfo[lastPatt] += ' ' + line

        if len(searchObj.groups() ) > 0:
          expInfo[patt] = searchObj.group(1)

      patt = self.expPattsList[-1]

      for line in expDetLines[lineNum:]:
        searchObj = self.emptyLinePatt.search(line)

        if searchObj:
          break
        else:
          expInfo[patt] += ' ' + line.strip()

      #for k in expInfo.keys():
      #  print 'DATA2: [' + k + '] [' + str(expInfo[k]) + ']'

    return expInfo


  def getSoftwareMethod(self, softwareName, methodName=None, authors=None, remarks=None):
  
    searchObj = self.softVersPatt.search(softwareName)

    if searchObj:
      softName = searchObj.group(1)
      versNum  = searchObj.group(2)

    else:
      softName = softwareName
      versNum  = 'any'

    software = self.ccpnProject.currentMethodStore.findFirstSoftware(name=softName)

    if not software:
      software = self.ccpnProject.currentMethodStore.findFirstSoftware(name=softwareName)

    if not software:
      software = self.ccpnProject.currentMethodStore.newSoftware(name=softName, version=versNum) # NEW

      if DEBUG:
        print 'NEW SOFTWARE: [%s]' % software

    method = None

    if methodName:
      method = software.findFirstMethod(name=methodName)

      if not method and methodName != 'NULL':
        method = self.ccpnProject.currentMethodStore.newMethod(name=methodName, software=software) # NEW

        if DEBUG:
          print 'NEW METHOD: [%s]' % method

    if authors:
      software.vendorName = authors

    if remarks:
      software.details = remarks

    return (software, method)


  def getVaguePerson(self,authorName, makePerson=True):

    authorName = authorName.strip()

    authorComps = authorName.split('.')
    if len(authorComps) == 1:
      authorComps = authorName.split()

    familyName = authorComps[-1]

    person = self.findVaguePerson(familyName,authorComps[0])

    if DEBUG:
      print 'FOUND PERSON: [%s]' % person

    if person:

      if len(person.givenName) == 1 and len(authorComps[0]) > 1:
        person.setGivenName(authorComps[0])

        if DEBUG:
          print 'SET PERSON NAME: [%s] [%s]' % (person, authorComps[0])

    elif makePerson:

      keywds = {'familyName': familyName.capitalize()}

      if len(authorComps) > 1:
        keywds['givenName'] = authorComps[0].capitalize()

      if len(authorComps) > 2:
        keywds['middleInitials'] = authorComps[1:-1]

      person = self.ccpnProject.currentAffiliationStore.newPerson(**keywds) # NEW

      if DEBUG:
        print 'NEW PERSON: [%s]' % person

    return person


  def findVaguePerson(self,familyName,firstNameGuess):

    person = None

    for tperson in self.ccpnProject.currentAffiliationStore.sortedPersons():
      if tperson.familyName and tperson and tperson.familyName.upper() == familyName.upper():
        if (tperson.givenName and firstNameGuess and \
            (tperson.givenName.upper() == firstNameGuess.upper() or \
             tperson.givenName[0].upper() == firstNameGuess[0].upper()) ) or \
             tperson.firstInitial == firstNameGuess:
          person = tperson
          break

    return person


  def getTopLevelStore(self, memopsRoot, className, name=None, defaultName='readPdb'):

    store = None

    if name is None:
      store = getattr(memopsRoot, 'current'+className) or \
              getattr(memopsRoot, 'findFirst'+className)()

      if not store:
        store = getattr(memopsRoot, 'new'+className)(name=defaultName)

        if DEBUG:
          print 'NEW STORE: [%s] [%s] [%s]' % (className, store, defaultName)

    else:
      store = getattr(memopsRoot, 'findFirst'+className)(name=name)

      if not store:
        store = getattr(memopsRoot, 'new'+className)(name=name)

        if DEBUG:
          print 'NEW STORE: [%s] [%s] [%s]' % (className, store, name)

    return store


  def getTopLevelNSStore(self, memopsRoot, className, namingSystem=None, defaultNamingSystem='readPdb'):

    store = None

    if namingSystem is None:
      store = getattr(memopsRoot, 'current'+className) or \
              getattr(memopsRoot, 'findFirst'+className)()

      if not store:
        store = getattr(memopsRoot, 'new'+className)(namingSystem=defaultNamingSystem)

        if DEBUG:
          print 'NEW STORE: [%s] [%s] [%s]' % (className, store, defaultNamingSystem)

    else:
      store = getattr(memopsRoot, 'findFirst'+className)(namingSystem=namingSystem)

      if not store:
        store = getattr(memopsRoot, 'new'+className)(namingSystem=namingSystem)

        if DEBUG:
          print 'NEW STORE: [%s] [%s] [%s]' % (className, store, namingSystem)

    return store


if __name__ == '__main__':

  # Add error checking here.

  pdbFileName = sys.argv[1]

  autodep = ReadPdb(pdbFileName)

  ccpnProject = autodep.ccpnProject

  ccpnProject.saveModified()
