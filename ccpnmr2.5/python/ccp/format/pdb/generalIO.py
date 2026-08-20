
"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for PDB files

Copyright (C) 2005-2006 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


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
import os

from ccp.format.general.formatIO import FormatFile
from ccp.format.general.Constants import defaultMolCode

# Import general functions
from memops.universal.Util import returnInt
from memops.universal.Io import getTopDirectory

from ccp.format.pdb.cifCodeRedirect import redirectDict

#####################
# Class definitions #
#####################

class PdbGenericFile(FormatFile):

  def setGeneric(self):
    
    self.format = 'pdb'
    self.defaultMolCode = defaultMolCode
    
    self.version = None

class PdbFile(PdbGenericFile):

  """
  
  PDB file parser. Written by Wim Vranken, additions by Jules Jacobson (JJ) from the EBI.
  
  """

  def initialize(self):

    self.chains = []
    self.molIdToChains = {}

    self.hetGroups = []
    self.modelCoordinates = {}
    
    self.code = None
    self.title = None

    self.ssBonds = []
    self.links = []
    self.dbRefs = []

    self.seqStrucInfo = {}
    
    self.headerVars = {}
    self.otherCitations = {}
    self.refineInfo = {}
    self.remarks = {}

    self.sourceInfo = []

    self.helixCodes = {
        1:  'Right-handed alpha',
        2:  'Right-handed omega',
        3:  'Right-handed pi',
        4:  'Right-handed gamma',
        5:  'Right-handed 310',
        6:  'Left-handed alpha',
        7:  'Left-handed omega',
        8:  'Left-handed gamma',
        9:  '27 ribbon/helix',
        10: 'Polyproline'
    }
    
    self.sourceVarInfo = {
      'SYNTHETIC':                  'synthetic',
      'FRAGMENT':                   'fragment',       # Domain or region of protein
      'ORGANISM_SCIENTIFIC':        'organismScientific',  
      'ORGANISM_COMMON':            'organismCommon',
      'STRAIN':                     'strain',
      'VARIANT':                    'variant',
      'CELL_LINE':                  'cellLine',
      'ATCC':                       'atcc',
      'ORGAN':                      'organ',
      'TISSUE':                     'tissue',
      'CELL':                       'cell',
      'ORGANELLE':                  'organelle',
      'GENE':                       'gene',
      'SECRETION':                  'secretion',
      'CELLULAR_LOCATION':          'cellularLocation',
      'ORGANISM_TAXID':             'taxId',
      'EXPRESSION_SYSTEM':          'expressionSystem',
      'EXPRESSION_SYSTEM_COMMON':   'expressionSystemCommon',
      'EXPRESSION_SYSTEM_STRAIN':   'expressionSystemStrain',
      'EXPRESSION_SYSTEM_VARIANT':  'expressionSystemVariant',
      'EXPRESSION_SYSTEM_CELL_LINE':'expressionSystemCellLine',
      'EXPRESSION_SYSTEM_ATCC_NUMBER':'expressionSystemAtccNumber',
      'EXPRESSION_SYSTEM_ORGAN':    'expressionSystemOrgan',
      'EXPRESSION_SYSTEM_TISSUE':   'expressionSystemTissue',
      'EXPRESSION_SYSTEM_CELL':     'expressionSystemCell',
      'EXPRESSION_SYSTEM_ORGANELLE':'expressionSystemOrganelle',
      'EXPRESSION_SYSTEM_CELLULAR_LOCATION':'expressionSystemCellularLocation',
      'EXPRESSION_SYSTEM_VECTOR_TYPE':'expressionSystemVectorType',
      'EXPRESSION_SYSTEM_VECTOR':   'expressionSystemVector',
      'EXPRESSION_SYSTEM_PLASMID':  'expressionSystemPlasmid',
      'EXPRESSION_SYSTEM_GENE':     'expressionSystemGene',
      'EXPRESSION_SYSTEM_TAXID':    'expressionSystemTaxId',
      'OTHER_DETAILS':              'expressionSystemOtherDetails'
      }
    
  def setCode(self,code):
  
    self.code = code
    
  def convertResName(self,resName):
  
    resName = resName.strip()
    
    if redirectDict.has_key(resName):
      resName = redirectDict[resName]
      
    return resName
    
  def setMolIdToChains(self,molId,chain):
    if not self.molIdToChains.has_key(molId):
      self.molIdToChains[molId] = []
    self.molIdToChains[molId].append(chain)

  def setSSBond(self,resName1,chainId1,seqCode1,iCode1,resName2,chainId2,seqCode2,iCode2):
      
    self.ssBonds.append(((resName1,chainId1,returnInt(seqCode1),iCode1),(resName2,chainId2,returnInt(seqCode2),iCode2)))

  def setLink(self,atomName1,resName1,chainId1,seqCode1,iCode1,atomName2,resName2,chainId2,seqCode2,iCode2):
      
    self.links.append(((atomName1,resName1,chainId1,returnInt(seqCode1),iCode1),(atomName2,resName2,chainId2,returnInt(seqCode2),iCode2)))

  def setDbRef(self,pdbCode,chainId,seqBegin,insertBegin,seqEnd,insertEnd,dataBase,dbAccession,dbIdCode,dbSeqBegin,idbnsBeg,dbSeqEnd,dbinsEnd):

    self.dbRefs.append( (pdbCode,chainId,seqBegin,insertBegin,seqEnd,insertEnd,dataBase,dbAccession,dbIdCode,dbSeqBegin,idbnsBeg,dbSeqEnd,dbinsEnd) )

  def readHeader(self,verbose = 0):

    #
    # Read header lines into memory
    #

    fin = open(self.name, 'rU')
    allLines = fin.readlines()
    fin.close()

    headerLines = []
    
    for line in allLines:

      if not self.patt[self.format + 'AllAtom'].search(line):
        headerLines.append(line)
      else:
        break

    #
    # Determine whether PDB file, also determine if in V2.0 or higher format
    #

    isPdb = 0
    isNew = 0
    chainNum = -1
    setChainDict = {}
    title = ""

    for line in headerLines:
      if self.patt[self.format + 'Header'].search(line):
        # PDB header present, get info
        isPdb = 1
        pdbCode = line[62:66]
        self.setCode(pdbCode)
      
      # Title added by JJ
      elif self.patt[self.format + 'Title'].search(line):
        tmpLine = line[10:]
        if tmpLine:
          tmpLine = tmpLine.strip()
        title += ' ' + tmpLine
        self.title = title

      elif self.patt[self.format + 'Remark4'].search(line):
        isNew = 1
        
        # Also now get version number for determining which atom naming
        # system to use
        
        versionSearch = self.patt[self.format + 'Remark4'].search(line)
        
        self.version = versionSearch.group(2)

    #
    # Start reading file according to format
    #  

    if isPdb == 1:

      # TODO: read in disulfide bonds from SSBOND line

      oldSerNum = -1
      synonym = None

      if isNew == 1:
        molId = None
        highMolId = -1
      else:
        molId = 1
        highMolId = molId

      molName = None
      tempMolName = ''

      lineNum = 0
      headerLineNum = len(headerLines)
      oldChainId = None

      # While is necessary to handle line jumps (when info for current line
      # is continued on next

      while lineNum < headerLineNum:

        line = headerLines[lineNum]

        if self.patt[self.format + 'Sequence'].search(line):
          # Polymer sequence
          
          if self.version in ('3.20',):
            offset = 1
            chainId = line[11+offset:12+offset]
            
            # Seem to be inconsistencies in files... hack to fix this
            if not chainId.strip():
              chainId = line[11:12]
            
          else:
            offset = 0
            chainId = line[11:12]

          serNum = returnInt(line[7:10+offset])
          numRes = line[13+offset:17+offset]                # Is repeated
          resNames = line[18+offset:70+offset].split()
 
          if serNum != oldSerNum + 1 or chainId != oldChainId:
            found = 0
            oldChainId = chainId

            # New chain, search which one
            for chain in self.chains:
              if chain.chainId == chainId:
                found = 1
                break

            # Exception to the chain stuff - NULL can be chain N (?? in 1arj)
            if not found and isNew and chainId == 'N':
              for chain in self.chains:
                if chain.chainId == self.defaultMolCode and not chain.numRes:
                  chain.chainId = chainId
                  chain.refChainId = chainId
                  found = 1
                  break
                  
            if not found:

              if isNew:
                # Error for new PDB files

                if not molName:
                  molName = tempMolName
                
                if chainId != self.defaultMolCode and len(self.chains) == 1 and self.chains[0].chainId == self.defaultMolCode and not self.chains[0].resNames:
                
                  # Use existing chain in case a NULL one is available...
                  print "Warning: replacing NULL chain by chain Id '%s' in pdb file %s." % (chainId,self.name)
                  
                  chain = self.chains[0]
                  chain.chainId = chainId
                  chain.refChainId = chainId
                  molId = chain.molId

                else:        
                  print "Error: unknown chain Id '%s' in pdb file %s. Created from scratch." % (chainId,self.name)
                  self.chains.append(PdbPolymer(self,chainId,molId,molName))
                  chain = self.chains[-1]
                  
                  self.setMolIdToChains(molId,chain)
                  
                  molId = molId + 1
                  highMolId = molId

              else:
                # For old PDB files use chain info from seqres
                # Molecule sequences are checked later (so same seq refers
                # to same molecule)
                if not molName:
                  molName = tempMolName
                
                self.chains.append(PdbPolymer(self,chainId,molId,molName))
                chain = self.chains[-1]
                
                self.setMolIdToChains(molId,chain)
                
                molId = molId + 1
                highMolId = molId

            chain.addNumRes(numRes)
            
          for tResName in resNames:
            chain.addResName(self.convertResName(tResName))

          oldSerNum = serNum

        elif self.patt[self.format + 'HetGroup'].search(line):
          # Hetero group (non standard stuff)

          hetId = self.convertResName(line[6:10])
          chainId = line[12:13]
          seqNum = line[13:17]
          iCode = line[17:18]
          numHetAtoms = line[20:25]
          description = line[30:70]
          
          # Need to reset chainId if listed as part of seqRes!
          if not chainId.strip():
            for chain in self.chains:
              if hetId in chain.resNames:
                print "  Warning: matching HET group to SEQRES chain ID..."
                chainId = chain.chainId

          self.hetGroups.append(PdbHetGroup(hetId,chainId,seqNum,iCode,numHetAtoms,description))

        elif self.patt[self.format + 'HetFormula'].search(line):
          # More het group info

          compNum = line[8:10]       # Component number (currently ignored). Can link hetgroup to seqres block
          hetId = self.convertResName(line[12:15])
          continuation = line[16:18]
          asterisk = line[16:18]     # Is * for water
          formula = line[18:70]

          found = 0
          
          hetGroupsFound = []

          for hetGroup in self.hetGroups:
            if hetGroup.Id == hetId:
              hetGroupsFound.append(hetGroup)
              found = 1

          if found == 0:
            print "Warning: PDB code %s, added hetero group '%s', formula %s from FORMUL line." % (self.code,hetId,formula)
            
            # in this case have to parse the HETATOM lines to find the correct chainId,seqNum,iCode,numHetAtoms
            hetIdsFound = []
            hetIdAtoms = {}
            
            for coordLine in allLines[headerLineNum:]:
              if self.patt[self.format + 'HetAtom'].search(coordLine):
                resName = self.convertResName(coordLine[17:20])

                if resName == hetId:
                  chainId = coordLine[21:22]
                  seqCode = coordLine[22:26]
                  insertionCode = coordLine[26:27]
                  
                  hetIdInfo = (chainId,seqCode,insertionCode)
                  
                  if hetIdInfo not in hetIdsFound:
                    hetIdsFound.append(hetIdInfo)
                    hetIdAtoms[hetIdInfo] = 0
                  
                  hetIdAtoms[hetIdInfo] += 1                

              if self.patt[self.format + 'NewModel'].search(line):
                break

            for (chainId,seqCode,insertionCode) in hetIdsFound:
              numAtoms = hetIdAtoms[(chainId,seqCode,insertionCode)]
              self.hetGroups.append(PdbHetGroup(hetId,chainId,seqCode,insertionCode,numAtoms,None))
              hetGroupsFound.append(self.hetGroups[-1])
            
          for hetGroup in hetGroupsFound:
            hetGroup.addFormula(formula)

        elif self.patt[self.format + 'SsBond'].search(line):
          # Disulfide bond info

          serNum = returnInt(line[7:10])

          resName1 = self.convertResName(line[11:14])
          chainId1 = line[15:16]
          seqCode1 = line[17:21].strip()
          iCode1 = line[21:22]

          resName2 = self.convertResName(line[25:28])
          chainId2 = line[29:30]
          seqCode2 = line[31:35].strip()
          iCode2 = line[35:36]

          # Set SS bonds between atoms after reading in atom coordinates.
          self.setSSBond(resName1,chainId1,seqCode1,iCode1,resName2,chainId2,seqCode2,iCode2)

        elif self.patt[self.format + 'Link'].search(line):
          # Linked residue info - useful to group chemComps into chain...
          # Lists atoms bonded between HET groups or between HET groups and normal residues

          atomName1 = line[12:16].strip()
          resName1 = self.convertResName(line[17:20])
          chainId1 = line[21:22]
          seqCode1 = returnInt(line[22:26])
          iCode1 = line[26:27]

          atomName2 = line[42:46].strip()
          resName2 = self.convertResName(line[47:50])
          chainId2 = line[51:52]
          seqCode2 = returnInt(line[52:56])
          iCode2 = line[56:57]

          # Set SS bonds between atoms after reading in atom coordinates.
          self.setLink(atomName1,resName1,chainId1,seqCode1,iCode1,atomName2,resName2,chainId2,seqCode2,iCode2)

        elif self.patt[self.format + 'DbReference'].search(line):
          # Cross-reference between PDB sequence and other database entries.
          # Not always present (only for peptide with 10 or more residues)
          # Currently ignored

          pdbCode = line[7:11]
          chainId = line[12]

          seqBegin = line[14:18]
          insertBegin = line[18]

          seqEnd = line[20:24]
          insertEnd = line[24]

          dataBase = line[26:32]
          dbAccession = line[33:41]
          dbIdCode = line[42:54]

          dbSeqBegin = line[55:60]
          idbnsBeg = line[60]

          dbSeqEnd = line[62:67]
          dbinsEnd = line[67]

          #print 'PDB: [' + pdbCode + '] [' + chainId + '] [' + seqBegin + '] [' + insertBegin + '] [' + seqEnd + '] [' + insertEnd + '] [' + dataBase + '] [' + dbAccession + '] [' + dbIdCode + '] [' + dbSeqBegin + '] [' + idbnsBeg + '] [' + dbSeqEnd + '] [' + dbinsEnd + ']'

          self.setDbRef(pdbCode,chainId,seqBegin,insertBegin,seqEnd,insertEnd,dataBase,dbAccession,dbIdCode,dbSeqBegin,idbnsBeg,dbSeqEnd,dbinsEnd)

        else:

          #
          # File level info
          #

          for searchInfo in ('Keywds','ExpData','Authors'):

            searchObj = self.patt[self.format + searchInfo].search(line[:70])

            if searchObj:

              otherText = searchObj.group(1)
              otherText = otherText.strip()

              if not self.headerVars.has_key(searchInfo):
                self.headerVars[searchInfo] = ""
              else:
                self.headerVars[searchInfo] += ' '

              self.headerVars[searchInfo] += otherText
              
          #
          # Journal
          #

          searchObj = self.patt[self.format + 'Journal'].search(line[:80])
          
          if searchObj:

            infoType = searchObj.group(1).strip()
            otherText = searchObj.group(2).strip()

            if not self.headerVars.has_key('Journal'):
              self.headerVars['Journal'] = {}
              
            self.handleJournalInfo(self.headerVars['Journal'],infoType,otherText)
                
          #
          # Remarks
          #

          searchObj = self.patt[self.format + 'Remarks'].search(line[:80])

          if searchObj:
            
            remarkNum = searchObj.group(1)
            otherText = searchObj.group(2).strip()

            if otherText:
              if remarkNum == '1':
                refIdSearch = self.patt[self.format + 'Reference'].search(otherText)
                if refIdSearch:
                  refId = int(refIdSearch.group(1).strip())
                  self.otherCitations[refId] = {}
                else:
                  searchObj = self.patt[self.format + 'ReferenceJournal'].search(otherText)
                  if searchObj:
                    infoType = searchObj.group(1).strip()
                    otherText = searchObj.group(2).strip()

                    self.handleJournalInfo(self.otherCitations[refId],infoType,otherText)

              else:
                
                if not self.remarks.has_key(remarkNum):
                  self.remarks[remarkNum] = ""

                if line[12:15] != '   ' and len(self.remarks[remarkNum]) > 0 and self.remarks[remarkNum][-1] != self.newline:
                  self.remarks[remarkNum] += self.newline

                self.remarks[remarkNum] += otherText

                if otherText[-1] not in ('(', ')'):
                  self.remarks[remarkNum] += self.newline

          searchObj = self.patt[self.format + 'Source'].search(line)
            
          if searchObj:
            
            otherText = searchObj.group(1).strip()

            searchTextMol = 'MOL_ID:'

            if otherText[:len(searchTextMol)] == searchTextMol:
              molIdText = otherText[len(searchTextMol):]
              molIdText = molIdText.replace(';','')
              molIdText = molIdText.strip()
              molId = returnInt(molIdText)

            searchTuple = tuple(self.sourceVarInfo.keys() )
            # was this for some reason:
            #searchTuple = ('ORGANISM_SCIENTIFIC','ORGANISM_COMMON','EXPRESSION_SYSTEM','EXPRESSION_SYSTEM_STRAIN') + tuple(self.sourceVarInfo.keys())

            if not hasattr(self, 'oldSearchText'):
              self.oldSearchText = ''

            for searchText in searchTuple:
              if otherText[:len(searchTextMol)] == searchTextMol:
                break

              if otherText[:len(searchText)] == searchText and otherText[len(searchText)] == ':':
                infoText = otherText[len(searchText)+2:]
                infoText = infoText.replace(';','')
                infoText = infoText.strip()

                if len(self.sourceInfo) == molId - 1:
                  self.sourceInfo.append({})

                self.sourceInfo[molId-1][searchText] = infoText
                self.oldSearchText = searchText
                break

            else:
              if self.oldSearchText:
                infoText = line[10:]
                infoText = infoText.replace(';','')
                infoText = infoText.strip()

                if len(self.sourceInfo) == molId - 1:
                  self.sourceInfo.append({})

                #print 'INFO: [%s] [%s]' % (self.sourceInfo[molId-1][self.oldSearchText], infoText)

                if self.sourceInfo[molId-1][self.oldSearchText][-1] in ('-', '.', '(', ')', '/', '\\', '%'):
                  self.sourceInfo[molId-1][self.oldSearchText] += '' + infoText
                else:
                  self.sourceInfo[molId-1][self.oldSearchText] += ' ' + infoText

        #
        # Information ONLY from new pdb files (V2.0 or higher)
        #

        if isNew == 1:
          
          # Some lines are multiple, sometimes there is
          # a ; at the end (new files), but not always (old ones).
          
          cmpSearchObj = self.patt[self.format + 'Compound'].search(line[:70])

          if cmpSearchObj:
            # Polymer info (stored in chains variable)
            line = line[:70].replace(cmpSearchObj.group(),'')
            colonSplit = line.strip().split(":") # Get rid of newline, get code
            code = colonSplit[0]
            value = ":".join(colonSplit[1:])
            value = value.strip()

            addLineNum = 1
            while (1):
              # Check if next line is part of this one
              addLine = headerLines[lineNum + addLineNum]

              if self.patt[self.format + 'Compound'].search(addLine):
                if self.patt['colon'].search(addLine):
                  lineNum += addLineNum - 1
                  break
                else:
                  #value += " " + string.strip(line[10:70])
                  value += " " + addLine[10:70].strip()
                  addLineNum += 1
              else:
                break

            value = value.replace(";",'')	         # Bit dangerous                

            if code == 'CHAIN':
              # Chains present for this molecule (comma-separated)
              chainCodes = value.split(',')
              chainNum = len(chainCodes)

              if chainNum == 1 and chainCodes[0] == 'NULL':
                # New molecule, one chain
                self.chains.append(PdbPolymer(self,self.defaultMolCode,molId,molName))
                self.setMolIdToChains(molId,self.chains[-1])

                if synonym:
                  self.chains[-1].addValue('SYNONYM',synonym)

                chainNum = 1

              else:
                for chainCode in chainCodes:
                  chainCode = chainCode.strip()
                  # New molecule, multiple chains
                  self.chains.append(PdbPolymer(self,chainCode,molId,molName))
                  self.setMolIdToChains(molId,self.chains[-1])

              if setChainDict:
              
                for code in setChainDict.keys():
                
                  for num in range(0,chainNum):
                    self.chains[-(num+1)].addValue(code,setChainDict[code])
                    
            elif code == 'MOL_ID':
              # Numbers molecule
              molId = returnInt(value)

              if molId > highMolId:
                highMolId = molId

            elif code == 'MOLECULE':
              # Molecule name
              molName = value

            elif len(self.chains) == 0 and code == 'SYNONYM':
              # Synonym (can be listed after CHAIN?)
              synonym = value

            else:
              if chainNum == -1:
                # Chain not done yet - keep track
                setChainDict[code] = value
                
              else:
                for num in range(0,chainNum):
                  self.chains[-(num+1)].addValue(code,value)

          elif self.patt[self.format + 'DbReference'].search(line):
            # Cross-reference between PDB sequence and other database entries.
            # Not always present (only for peptide with 10 or more residues)
            # Currently ignored

            pass

          elif self.patt[self.format + 'SequenceChange'].search(line):

            # Conflicts between db reference and actual sequence of molecule
            # Currently ignored

            pass

          elif self.patt[self.format + 'HetName'].search(line):  # This also in old PDB?
            # More het group info

            continuation = line[8:10]
            hetId = self.convertResName(line[10:14])
            chemicalName = line[15:70]

            for hetGroup in self.hetGroups:
              if hetGroup.Id == hetId:
                hetGroup.addChemicalName(chemicalName)


          elif self.patt[self.format + 'HetSynonym'].search(line):  # This also in old PDB?
            # More het group info

            continuation = line[8:10]
            hetId = self.convertResName(line[10:14])
            synonym = line[15:70]

            for hetGroup in self.hetGroups:
              if hetGroup.Id == hetId:
                hetGroup.addSynonym(synonym)
          
          
          elif self.patt[self.format + 'SecStruc'].search(line):
            
            secStrucType = line[:6].strip()
            
            serial = returnInt(line[7:10])
            ssId = line[11:14].strip()
            
            specificInfo = {}
            
            if secStrucType == 'HELIX':
            
              initResName = line[15:18].strip()
              initChainID = line[19:20]
              initSeqNum = returnInt(line[21:25])
              initInsertCode = line[25:26]
              
              endResName = line[27:30].strip()
              endChainID = line[31:32]
              endSeqNum = returnInt(line[33:37])
              endInsertCode = line[37:38]

              helixClassId = returnInt(line[38:40])
              if self.helixCodes.has_key(helixClassId):
                specificInfo['class'] = self.helixCodes[helixClassId]
              else:
                print "  Warning: no helix class found for identifier '%s'!" % str(helixClassId)
                specificInfo['class'] = None
              specificInfo['comment'] = line[40:70].strip()
              specificInfo['length'] = returnInt(line[71:76])
              
            elif secStrucType == 'SHEET':

              numStrands = returnInt(line[14:16])
              
              initResName = line[17:20].strip()
              initChainID = line[21:22]
              initSeqNum = returnInt(line[22:26])
              initInsertCode = line[26:27]
              endResName = line[28:31].strip()
              endChainID = line[32:33]
              endSeqNum = returnInt(line[33:37])
              endInsertCode = line[37:38]
              
              specificInfo['sense'] = returnInt(line[38:40]) # 0 if first, 1 if parallel, -1 if anti-parallel with respect to previous
              
              """
              TODO: Registration information on hydrogen bonding...
              
              42 - 45   Atom            curAtom         Registration. Atom name in current
                                                        strand.
              46 - 48   Residue name    curResName      Registration. Residue name in
                                                        current strand.
              50        Character       curChainId      Registration. Chain identifier in
                                                        current strand.
              51 - 54   Integer         curResSeq       Registration. Residue sequence
                                                        number in current strand.
              55        AChar           curICode        Registration. Insertion code in
                                                        current strand.
              57 - 60   Atom            prevAtom        Registration. Atom name in
                                                        previous strand.
              61 - 63   Residue name    prevResName     Registration. Residue name in
                                                        previous strand.
              65        Character       prevChainId     Registration. Chain identifier in
                                                        previous strand.
              66 - 69   Integer         prevResSeq      Registration. Residue sequence
                                                        number in previous strand.
              70        AChar           prevICode       Registration. Insertion code in
                                                        previous strand.
              """

            elif secStrucType == 'TURN':

              initResName = line[15:18].strip()
              initChainID = line[19:20]
              initSeqNum = returnInt(line[20:24])
              initInsertCode = line[24:25]
              
              endResName = line[26:29].strip()
              endChainID = line[30:31]
              endSeqNum = returnInt(line[31:35])
              endInsertCode = line[35:36]
              
              specificInfo['comment'] = line[40:70].strip()
              
              
            else:
              print "  Warning: unknown secondary structure type %s! Ignored" % secStrucType
            
            #
            # Now set the information for later reference...
            #
            
            if not self.seqStrucInfo.has_key(secStrucType):
              self.seqStrucInfo[secStrucType] = []
            
            self.seqStrucInfo[secStrucType].append((serial,ssId,
                                                   (initResName,initChainID,initSeqNum,initInsertCode),
                                                   (endResName,endChainID,endSeqNum,endInsertCode),
                                                   specificInfo))
            
          else:
            
            #
            # Other molecule level info
            #

            searchObj = self.patt[self.format + 'Source'].search(line)
            
            if searchObj:
            
              otherText = searchObj.group(1).strip()

              searchTextMol = 'MOL_ID:'
              if otherText[:len(searchTextMol)] == searchTextMol:
                molIdText = otherText[len(searchTextMol):]
                molIdText = molIdText.replace(';','')
                molIdText = molIdText.strip()
                molId = returnInt(molIdText)

              searchTuple = tuple(self.sourceVarInfo.keys() )
              #searchTuple = ('ORGANISM_SCIENTIFIC','ORGANISM_COMMON','EXPRESSION_SYSTEM','EXPRESSION_SYSTEM_STRAIN') + tuple(self.sourceVarInfo.keys())

              if not self.molIdToChains.has_key(molId):
                print "  Error: no molIdToChains defined for MOL_ID %s" % str(molId)

              else:
                if not hasattr(self, 'oldSearchText'):
                  self.oldSearchText = ''
                for searchText in searchTuple:
                  if otherText[:len(searchTextMol)] == searchTextMol:
                    break

                  if otherText[:len(searchText)] == searchText and otherText[len(searchText)] == ':':
                    infoText = otherText[len(searchText)+2:]
                    infoText = infoText.replace(';','')
                    infoText = infoText.strip()
                    for chain in self.molIdToChains[molId]:
                      chain.headerVars[searchText] = infoText
                      chain.addSource(searchText,infoText) #JJ            

                    self.oldSearchText = searchText
                    break

                else:
                  infoText = line[10:]
                  infoText = infoText.replace(';','')
                  infoText = infoText.strip()

                  for chain in self.molIdToChains[molId]:
                    if self.oldSearchText:
                      if chain.headerVars.has_key(self.oldSearchText):
                        chain.headerVars[self.oldSearchText] += ' ' + infoText
                        chain.addSourceExtra(self.oldSearchText,infoText)
                      else:
                        chain.headerVars[self.oldSearchText] = infoText
                        chain.addSource(self.oldSearchText,infoText)

        #
        # Information ONLY from old pdb files 
        #

        else:

          cmpSearchObj = self.patt[self.format + 'Compound'].search(line[:70])

          if cmpSearchObj:
            # Polymer info (stored in chains variable)
            line = line[:70].replace(cmpSearchObj.group(),'')
            
            tempMolName += line.strip()

            if line.count('MOLECULE'):
              colonSplit = line.strip().split(":") # Get rid of newline, get code
              code = colonSplit[0]
              value = ":".join(colonSplit[1:])
              value = value.strip()
              if value[-1] == ';':
                value = value[:-1]
              
              molName = value

        lineNum += 1

      if isNew == 0:
        # Compare sequences for chains read in from SEQRES (old PDB)
        # Reset molId if same molecules, add info to molName if multiple molecules

        multipleMol = 0

        for i in range(0,len(self.chains) - 1):
          chain = self.chains[i]
          for j in range(i + 1,len(self.chains)):
            otherChain = self.chains[j]

            if chain.resNames == otherChain.resNames:
              otherChain.molId = chain.molId
            else:
              multipleMol = 1

        # Add (mol <n>) to molName (should be unique!)
        if multipleMol == 1:
          for i in range(0,len(self.chains)):
            self.chains[i].molName = self.chains[i].molName + (" (mol %d)" % (i+1))

      else:
        #
        # For new molecules still possible that one molecule name with multiple
        # chains with different sequence (eg insulin 1a7f)
        # Make sure molName and molId are different in this case
        #

        for i in range(0,len(self.chains) - 1):

          chain = self.chains[i]
          chainChanged = 0

          for j in range(i + 1,len(self.chains)):
            otherChain = self.chains[j]

            if chain.molName == otherChain.molName and chain.molId == otherChain.molId:
              if chain.resNames != otherChain.resNames:
                otherChain.molName = otherChain.molName + " (chain %s)" % otherChain.chainId
                highMolId += 1
                otherChain.molId = highMolId

                if chainChanged == 0:
                  chain.molName = chain.molName  + " (chain %s)" % chain.chainId


    else:
      print "Error: %s is not a file originating from PDB. Cannot be read." % self.name

#    for key in self.remarks.keys():
#      self.remarks[key] += self.newline
#      print 'KEY: [%s]' % self.remarks[key]

            
  def handleJournalInfo(self,infoDict,infoType,otherText):

    otherText = otherText.strip()

    if not infoDict.has_key(infoType):
      infoDict[infoType] = ""

    if infoType in ('TITL','AUTH','PUBL','EDIT'):
      if otherText[-1] == '-' or otherText[-1] == ':':
        infoDict[infoType] += otherText
      else:
        infoDict[infoType] += otherText + " "
    elif infoType in ('REF'):
      if not infoDict[infoType]:
        infoDict[infoType] = otherText
      else:
        infoDict[infoType] = infoDict[infoType][:-18].strip() + " " + otherText + infoDict[infoType][-18:]

class PdbPolymer:

  def __init__(self,parent,chainId,molId,molName):
  
    self.parent = parent

    self.chainId = chainId
    self.refChainId = chainId
    self.molId = molId
    self.molName = molName
    self.synonyms = []
    self.ecNums = []

    self.numRes = None
    self.resNames = []
    
    self.headerVars = {}
    
    for attrName in ('engineered','mutation','biologicalUnit','details') + tuple(self.parent.sourceVarInfo.values()):
      setattr(self,attrName,None)

  def addResName(self,resName):
  
    self.resNames.append(resName)
    
  def addNumRes(self,numRes):
  
    self.numRes = returnInt(numRes)

  def addValue(self,code,value):

    if code == 'SYNONYM':
      # Synonyms for this molecule (comma-separated)
      synonyms = value.split(',')
      for synonym in synonyms:
        self.synonyms.append(synonym.strip())

    elif code == 'EC':
      # Enzyme commission number (comma-separated)
      ecNums = value.split(',')
      for ecNum in ecNums:
        self.ecNums.append(ecNum.strip())

    elif code == 'ENGINEERED':
      # Made using recombinant technology or by chemical synthesis
      self.engineered = value
      
    elif code == 'MUTATION':
      # Describes mutations from wild type sequence
      self.mutation = value
      
    elif code == 'FRAGMENT':
      # Describes fragments from ???
      self.fragment = value
      
    elif code == 'BIOLOGICAL_UNIT':
      # If part of larger unit, describes whole functional unit
      self.biologicalUnit = value
      
    elif code == 'OTHER_DETAILS':
      # Additional comments.
      self.details = value

  def addSource(self,code,value):
  
    attrName = self.parent.sourceVarInfo[code]

    setattr(self,attrName,value)

  def addSourceExtra(self,code,value):
  
    attrName = self.parent.sourceVarInfo[code]

    sourceInfo = getattr(self,attrName)

    setattr(self,attrName,sourceInfo + ' ' + value)

class PdbHetGroup:

  def __init__(self,hetId,chainId,seqNum,iCode,numHetAtoms,description):
  
    self.Id = hetId
    self.chainId = chainId
    self.refChainId = chainId
    self.seqNum = returnInt(seqNum)
    self.iCode = iCode
    self.numHetAtoms = returnInt(numHetAtoms)
    
    if description:
      self.description = description.strip()
    else:
      self.description = None
    
    # Special case
    if hetId == 'HOH':
      self.chemicalName = 'water'
    else:
      self.chemicalName = ""
      
    self.synonyms = []
    self.formula = ""
    
  def addChemicalName(self,chemicalName):
    self.chemicalName += chemicalName.strip()
    
  def addSynonym(self,synonym):
    # Seems to be ; separated?
    self.synonyms.append(synonym.strip())
  
  def addFormula(self,formula):
    self.formula += formula.strip()

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = [
           # New style
           '../../reference/pdb/1hbb.pdb',
           '../../reference/pdb/1aa9.pdb',
           '../../reference/pdb/1b4c.pdb',
           '../../reference/pdb/1a7f.pdb',
           '../../reference/pdb/1a57.pdb',
           '../../reference/ccpNmr/jurgenBmrb/1klr/1klr.pdb',
           
           # This one is old style, no heteroatoms
           '../../rawdata/4/4/1/info.general/1shg.pdb',
           
           # This one is old style, has heteroatoms
           '../../reference/pdb/3xis.pdb',
           
           # This one is old style, has different chains
           '../../reference/pdb/4sbv.pdb',

           # Problems with sequence reading 
           '../../recoord/allDataTest/1ai0/1ai0.pdb',

           # Problems with sequence reading 
           '../../recoord/allDataTest/1ctl/1ctl.pdb',
           
           '../../reference/pdb/pdbnew/pdb-set-2/2fn2.pdb'
          ]

  
  for inFile in files[-1:]:
    
    pdbFile = PdbFile(os.path.join(getTopDirectory(),inFile))

    pdbFile.readHeader(verbose = 1)
    
    print pdbFile.headerVars
    print pdbFile.remarks
    
    print pdbFile.chains
    
    for chain in pdbFile.chains:
      print chain.chainId, chain.molId, chain.headerVars
      
    print pdbFile.hetGroups

    for hetGroup in pdbFile.hetGroups:
      alist = dir(hetGroup)
      for attr in alist:
        if attr in ['__init__','__doc__','__module__','addChemicalName','addFormula','addSynonym']:
          continue
        print attr, getattr(hetGroup,attr)
    
