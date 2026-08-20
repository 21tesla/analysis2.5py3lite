"""
======================COPYRIGHT/LICENSE START==========================

Constants.py: Useful constants for scripts in ccpnmr/format and subdirectories

Copyright (C) 2008-2010 Wim Vranken (European Bioinformatics Institute)

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

from ccp.format.general.Constants import defaultMolCode
from ccp.general.Constants import code1LetterToCcpCodeDict

import os

########################
# Variable definitions #
########################

coordinateImportSourceName = 'ccpnImport'

tagSep = '.'
newline = os.linesep
stereoText = " (stereo)"

#
# For constraints
#

distanceConstraintDefaultLowerLimit = 1.80

#
# For chemical shifts: if a value already exists for a given
# chemical shift list and trying to set a new one, then the the
# new value is ignored in the scripts (can still force setting
# this)
#

chemShiftAllowedDiff = 0.01

#######################
# Keyword definitions #
#######################

ccpNmr_kw = 'ccpNmr'
isAriaInput_kw = 'isAriaInput'
isAriaOutput_kw = 'isAriaOutput'
originalFormat_kw = 'originalFormat'
peakNum_kw = 'peakNum'
assign_kw = 'assign'
volume_kw = 'volume'
height_kw = 'height'
ccpNmrMapping_kw = 'ccpNmrMapping'
importGeneral_kw = 'importGeneral'
newResonances_kw = 'newResonances'
atomSerial_kw = 'atomSerial'
dateTimeStamp_kw = 'dateTimeStamp'
ioCycleTag_kw = 'ioCycleTag'
seqAssignTag_kw = 'seqAssign'
allResidueAtoms_kw = 'ALL_ATOMS'

#####################################################
# Definition of some frequency conversion constants #
#####################################################

from ccp.general.Constants import chemShiftRefRatios
INDIRECT_FREQ_CONV = {}
for tag in chemShiftRefRatios:
  INDIRECT_FREQ_CONV[tag] = 100.*chemShiftRefRatios[tag]
#INDIRECT_FREQ_CONV = {'1H': 100,
#                      '2H': 15.3506088,
#                      '13C': 25.1449530,
#                      '15N': 10.1329118,
#                      '31P': 40.4808636}

commonNmrActiveIsotopes = {'H': 1,
                           'C': 13,
                           'N': 15,
                           'P': 31}
                        

##########################################################
# Conversion of raw format attribute names to CCPN names #
##########################################################

rawToCcpnAttrNameDict = { 'point':    'position',
                          'pointDev': 'positionError',
                          'ppm':      'value',
                          'ppmDev':   'valueError'}

########################################
# Definition of chemical shift classes #
########################################

chemShiftClasses = {

  'protein':  {

      'N':  {

          '15N backbone':  [['N'], []],
          '15N sidechain': [[],    ['N']]

              },

      'C':  {

          '13C CO':        [['C'], []],
          '13C CA':        [['CA'],[]],
          '13C CB':        [['CB'],[]],
          '13C sidechain': [[],    ['C','CA','CB']]

              },

      'H':  {

          '1H HN':              [['H'], []],
          '1H HA':              [['HA','HA2','HA3'],[]],
          '1H sidechain (all)': [[],    ['H','HA','HA2','HA3']]

              }

      }
  }  

######################################################################################
# Definitions of all format lowercase (formatName) and uppercase (formatLabel) names #
######################################################################################

allFormatsDict = {'amber': 'Amber',
                  'ansig': 'Ansig',
                  'aqua': 'Aqua',
                  'aria': 'Aria',
                  'auremol': 'Auremol',
                  'autoAssign': 'AutoAssign',
                  'azara': 'Azara',
                  'bruker': 'Bruker',
                  'ccpNmr': 'CcpNmr',
                  'charmm': 'Charmm',
                  'cns': 'Cns',
                  'csi': 'Csi',
                  'concoord': 'Concoord',
                  'cosmos': 'Cosmos',
                  'cyana': 'Cyana',
                  'discover': 'Discover',
                  'dyana': 'Dyana',
                  'fasta': 'Fasta',
                  'felix': 'Felix',
                  'generic': 'Generic',
                  'gromacs': 'Gromacs',
                  'gromos': 'Gromos',
                  'mardigras': 'Mardigras',
                  'mars': 'Mars',
                  'marvin': 'Marvin',
                  'mmCif': 'MMCIF',
                  'module': 'Module',
                  'mol': 'Mol',
                  'mol2': 'Mol2',
                  'molmol': 'Molmol',
                  'monte': 'Monte',
                  'nmrDraw': 'NmrDraw',
                  'nmrPipe': 'NmrPipe',
                  'nmrStar': 'NmrStar',
                  'nmrView': 'NmrView',
                  'pales': 'Pales',
                  'pdb': 'Pdb',
                  'pipp': 'Pipp',
                  'pistachio': 'Pistachio',
                  'pronto': 'Pronto',
                  'pseudoPdb': 'PseudoPdb',
                  'regine': 'Regine',
                  'shiftx': 'Shiftx',
                  'shifty': 'Shifty',
                  'sparky': 'Sparky',
                  'talos': 'Talos',
                  'targetedAcquisition':'TargetedAcquisition',
                  'varian': 'Varian',
                  'xeasy': 'XEasy'}

#           
# Default naming systems for formats.
#
# Key is tuple with (formatName,formatVersion), value is naming system name (these are always uppercase)
#
                  
formatNamingDict = {('amber',None):      'PDB',
                    ('ansig',None):      'XPLOR',
                    ('aqua',None):       'XPLOR',  # TODO have to be able to set this!!
                    ('auremol',None):    'IUPAC',     # TODO: is this correct?
                    ('autoAssign',None): 'IUPAC',  # Is not relevant
                    ('aria',None):       'IUPAC',        # IS this correct?
                    ('ccpNmr',None):     'IUPAC',
                    ('charmm',None):     'XPLOR',
                    ('cns',None):        'XPLOR',
                    ('csi',None):        'IUPAC',
                    ('concoord',None):   'XPLOR',
                    ('cosmos',None):     'IUPAC',
                    ('cyana',None):      'CYANA2.1',   # This is the default naming system...
                    ('cyana','2.0'):     'DIANA',
                    ('cyana','2.1'):     'CYANA2.1',
                    ('diana',None):      'DIANA',
                    ('discover',None):   'MSI',
                    ('dyana',None):      'DIANA',
                    ('gromacs',None):    'IUPAC',
                    ('gromos',None):     'PDB',
                    ('mardigras',None):  'PDB',    # Probably not correct.
                    ('mars',None):       'IUPAC',  # Not really relevant - uses assignNames info for cs export
                    ('marvin',None):     'XPLOR',  # Runs in XPLOR, atom names come from there.
                    ('mmCif',None):      'PDB_REMED',
                    ('module',None):     'DIANA',
                    ('mol2',None):       'IUPAC',
                    ('molmol',None):     'MOLMOL',
                    ('nmrDraw',None):    'IUPAC',# TODO: is this correct??!
                    
                    ('nmrStar',None):    'AQUA',      # TODO: this might not be correct for non-constraint lists...
                    ('nmrStar','3.0'):   'PDB_REMED', # TODO: this might not be correct for non-constraint lists...
                    ('nmrStar','3.1'):   'PDB_REMED',
                    
                    ('nmrView',None):    'XPLOR',
                    ('pales',None):      'IUPAC',   # Pales is naming system agnostic TODO add other naming systems?
                    ('pdb',None):        'PDB',     # This is the default naming system
                    ('pdb','2.5'):       'IUPAC',   # This is an intermediate format...
                    ('pdb','3.0'):       'PDB_REMED',  
                    ('pdb','3.1'):       'PDB_REMED',  
                    ('pdb','3.2'):       'PDB_REMED',  
                                                
                    ('pipp',None):       'XPLOR',
                    ('pseudoPdb',None):  'PDB',
                    ('shiftx',None):     'IUPAC',
                    ('shifty',None):     'IUPAC',
                    ('sparky',None):     'DIANA',
                    ('talos',None):      'XPLOR', 
                    ('targetedAcquisition', None):'IUPAC',
                    ('xeasy',None):      'DIANA'}

# TODO: to get a full list have to do some more work... but this will do for now
namingSystemNames = []
for namingSystemName in formatNamingDict.values():
  if not namingSystemName in namingSystemNames:
    namingSystemNames.append(namingSystemName)
namingSystemNames.sort()
namingSystemNames.insert(0,'None')

formatComponents = {

  'single': [

     ['Sequence','sequence'],
     ['Chemical compounds','chemComps'],
     ['Coordinates','coordinates'],
     ['Proton exchange protection','hExchProtection'],
     ['Proton exchange rates','hExchRates'],
     ['Chemical shift anisotropy','shiftAnisotropy'],
     ['Chemical shifts','shifts'],
     ['T1 relaxation','t1Relax'],
     ['T1 rho relaxation','t1RhoRelax'],
     ['T2 relaxation','t2Relax'],
     ['Hetero NOEs','noe'],
     ['J couplings','jCoupling'],
     ['RDC values','rdc'],
     ['Peaks','peaks'],
     ['Peak assignments','peakAssignments'],
     ['Dihedral constraints','dihedralConstraints'],
     ['Distance constraints','distanceConstraints'],
     ['H-bond constraints','hBondConstraints'],
     ['Rdc constraints','rdcConstraints'],
     ['CSA constraints','csaConstraints'],
     ['J coupling constraints','jCouplingConstraints'],
     ['Acquisition parameters','acqPars'],
     ['Processing parameters','processingPars'],
     ['Processing scripts','processingScript'],
     ['People and citations','peopleAndCitations'],
     ['Assignment','assigment'],
     ['Project file','project']
     
            ],
            
  'combined': [
  
     ['Chemical shifts and peaks','shiftsAndPeaks']
     
              ]

}

# Used in export windows to determine what 'chainCode' info (e.g.
# segID for CNS) is necessary.

exportChainCodeDict = {'ansig': 'molecule name',
                       'auremol': 'molecule label',
                       'charmm':   'segID',
                       'cns':   'segID',
                       'concoord':   'segid',
                       'fasta': 'molecule name',
                       'nmrStar': 'Entity assembly (chain)',
                       'pdb': 'Chain ID'}


#
# Dictionary for project writing - which objects are required?
#

projectExportRequirements = {
                             
    'chain':     ('autoAssign','mars','talos'),
    'molSystem': ('sparky'),
    'shiftList': ('mars','talos'),
    'entry':     ('nmrStar','sparky')
    
}

############################
# Input/output definitions #
############################

"""
Following adapted from notes by Gary Thompson (February 2008):

 IoSetupList = [ IoGuiSetupListComponent1, IoGuiSetupListComponent2, ....]

    IoGuiSetupList = [
                      componentName,
                      componentReaderName,   # As used in DataFormat.py. Use '' for a lack of this operation).
                      componentWriterName,   # As used in DataFormat.py. Use '' for a lack of this operation).
                      guiSetupDict,          # A dictionary of guiSetups (see below) keyed by format-name.
                      inputGuiSetupLists,    # A list of guiSetupList (see below) defining default gui elements for import.
                      outputGuiSetupLists,   # A list of guiSetupLists (see below) defining default gui elements for export.
                      preExportDefinitionString, # String defining a function to call before exports are carried out.
                      postImportDefinitionString # String defining a function to call after imports are carried out.
                     ]
                     
guiSetup =  [
             default-file-name-string,
             importSupported,             # Is a boolean
             exportSupported,             # Is a boolean
             importGuiSetupList,          # This is format specific gui elements for import
             exportGuiSetupList           # This is format specific gui elements for export
            ]

guiSetupList =  [
                 guiElementTypeString,            # Widget name (SelectionList, FileButton, StringEntry, CheckButton, ...)
                 guiElementFcKeyword,             # Keyword used for DataFormat method calls (the value returned by the widget
                                                  # will be forwarded to the FC using this keyword).
                 guiElementLabel,                 # Text that appears in the GUI
                 guiElementDefaultValues,         # Default values for the widget
                 guiElementSpecialisationList     # Keywords to pass into the widget
                ]
"""

IoSetupList = [

   ['sequence',
    'readSequence',
    'writeSequence',
    {'amber':      ['amber.pdb',1,0,           [],[]],
     'ansig':      ['ansig.seq',1,1,           [],[]],
     'auremol':    ['auremol.seq',1,0,         [],[]],
     'autoAssign': ['table.aat',1,0,           [],[]],
     'charmm':     ['charmm.crd',1,0,          [],[]],
     'cns':        ['cnsOrXplor.pdb',1,0,      [],[]],
     'csi':        ['data.csi',1,0,            [],[]],
     'dyana':      ['dyana.seq',1,1,           [],[]],
     'cyana':      ['cyana.seq',1,1,           [['SelectionList','version','Cyana version:',['2.0','2.1'],None]],
                                               [['SelectionList','version','Cyana version:',['2.0','2.1'],None]]],
     'fasta':      ['fasta.seq',1,1,           [],[]],
     'gromos':     ['gromos.gro',1,0,          [],[]],
     'mars':       ['mars.seq',1,1,            [],[]],
     'mmCif':      ['project.cif',1,0,         [],[]],
     'mol2':       ['compound.mol2',1,0,       [],[]],
     'molmol':     ['molmol.pdb',1,0,          [],[]],
     'monte':      ['monte.seq',1,1,           [],[]],
     'nmrStar':    ['nmrStar.str',1,0,         [],[]],
     'nmrView':    ['nmrView.seq',1,1,         [],[]],
     'pdb':        ['code.pdb',1,0,            [],[]],
     'pistachio':  ['molecule.csa',1,0,        [],[]],
     'pronto':     ['ss.report',1,0,           [],[]],
     'pseudoPdb':  ['code.pdb',1,0,            [],[]],
     'shiftx':     ['shifts.tab',1,0,          [],[]],
     'sparky':     ['sparky.seq',1,1,          [],[]],
     'talos':      ['project.tab',1,0,         [],[]],
     'xeasy':      ['xeasy.seq',1,1,           [],[]]},
     [['FileButton','fileName','Sequence file to import:','Select file',None],
      ['SelectionListWithNone','molSystem','Use existing molecular system:','project.sortedMolSystems()',None],
      ['StringEntry','molName','Molecule name:','',{'width': 20}],
      ['IntEntry','numIdenticalChains','Number identical chains:',0,None],
      ['IntEntry','seqCodeStart','Sequence code start:',1,None],
      ['CheckButton','autoCreateChemComps','Automatically create chemComps if atom names available:',0,None],
      ['CheckButton','ignoreUnknownChemComps','Ignore unknown chemComps:',0,None],
      ['CheckButton','ignoreSingleOther','Ignore single residue chains:',0,None],
      ['CheckButton','autoLink','Fully automatic linking:',0,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],

    #'chainCodes':             (None,False,"List of chainCodes to be applied (default 'A','B','C',...)."),
    #'addNameMappings':        (None,False,'Dictionary that can be passed in to map atom names to CCPN atoms.'),
    #'ignoreResNames':         (None,False,"List of residue names that should be ignored when reading sequence (e.g. ['HOH'])."), 
    #'chemCompMappings':       (None,False,'Dictionary that can be passed in to map chemComp names from the format to CCPN chemComp names.'),
    #'chemCompPaths':          (None,False,'Dictionary with directories that contain chemComp information. Valid keys are chemCompData and chemCompDataLocal.'),

     [['MultiSelectionList','chains','Select chains to export:','project.sortedMolSystems().sortedChains()',None],
      ['FileButton','fileName','Name of export sequence file:','Select export file',None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',1,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['chemComps',
    'readChemComps',
    'writeChemComp',
    {'mol':     ['chemComp.mol',1,0,            [],[]],
     'mol2':    ['chemComp.mol2',1,1,           [],[]],
     'pdb':     ['code.pdb',1,0,                [],[]]},
     [['FileButton','fileName','File with chemical compound definition to import:','Select file',None],
      ['CheckButton','saveChemComp','Immediately save chemComp:',0,None],
      ['CheckButton','createAtomSets','Create atom sets:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],

    #'molTypes':                 (None,False,'Molecule types for the new chemComps'),
    #'ccpCodes':                 (None,False,'Unique CCP codes for the new chemComps'),
    #'chemCompPaths':            (None,False,'Dictionary with directories that contain chemComp information. Valid keys are chemCompData and chemCompDataLocal.'),

     [['SelectionList','chemCompVar','Select chemCompVar to export:','project.sortedMolSystems().sortedChains().sortedResidues().chemCompVar',None],
      ['FileButton','fileName','Name of export chemComp file:','Select export file',None],
      ['SelectionList','coordSystem','Default coordinate system:',['ideal','pdb'],None],
      ['CheckButton','writeLinkAtoms','Write out link atoms:',1,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['coordinates',
    'readCoordinates',
    'writeCoordinates',
    {'amber':     ['amber.pdb',1,0,           [],[]],
     'auremol':   ['auremol.pdb',1,1,         [],[]],
     'charmm':    ['charmm.crd',1,1,          [],[]],
     'cns':       ['cnsOrXplor.pdb',1,1,      [],[]],
     'cosmos':    ['cosmos.cod',0,1,          [],[]],
     'dyana':     ['dyana.cor',1,1,           [],[]],
     'cyana':     ['cyana.cor',1,1,           [['SelectionList','version','Cyana version:',['2.0','2.1'],None]],
                                              [['SelectionList','version','Cyana version:',['2.0','2.1'],None]]],
     'gromos':    ['gromos.gro',1,0,          [],[]],
     'mmCif':     ['project.cif',1,0,         [],[]],
     'molmol':    ['molmol.pdb',1,1,          [],[]],
     'pdb':       ['code.pdb',1,1,            [],[]],
     'pseudoPdb': ['code.pdb',1,1,            [['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
                                               ['CheckButton','swapFirstNumberAtom','Swap initial number in atom names:',0,None]],
                                              [['SelectionList','forceNamingSystemName','Use naming system:',namingSystemNames[1:],None]]]},
     [['FileMultiButton','fileNames','Coordinate files to import:','Select files',None],
      ['SelectionListWithNone','molSystem','Use existing molecular system:','project.sortedMolSystems()',None],
      ['SelectionListWithNone','structureEnsemble','Use existing structure ensemble:','project.currentNmrProject.sortedStructureEnsembles()',None],
      ['SelectionListWithNone','strucGen','Use existing structure generation:','project.currentNmrProject.sortedStructureGenerations()',None],
      ['IntEntry','maxNum','Maximum number of models to read:',9999,None],
      ['CheckButton','linkAtoms','Link unrecognized atoms:',1,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','matchToAtomNames','Match to internal atom names:',0,None],
      ['CheckButton','forceDefaultChainMapping','Use default chain mapping:',0,None],
      ['CheckButton','forceReadSequence','Force sequence read:',0,None],
      ['CheckButton','autoCreateChemComps','Automatically create chemComps if atom names available:',0,None],
      ['CheckButton','ignoreUnknownChemComps','Ignore unknown chemComps:',0,None],
      ['CheckButton','ignoreSingleOther','Ignore single residue molecules:',0,None],
      ['CheckButton','trackAtomSerials','Track original atom serials:',0,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],

    #'forceChainMappings':       (None,False,'Dictionary that can be passed in to automate the coordinate\n to data model chain/residue mapping.'),
    #'addNameMappings':          (None,False,'Dictionary that can be passed in to map coordinate atom names to CCPN atoms.'),
    #'ignoreResNames':           (None,False,"List of residue names that should be ignored when reading coordinates (e.g. ['HOH'])."), 
    #'chemCompPaths':            (None,False,'Dictionary with directories that contain chemComp information. Valid keys are chemCompData and chemCompDataLocal.'),
     
     [['MultiSelectionList','structures','Select structures to export:','project.sortedMolSystems().sortedStructureEnsembles().sortedModels()',None],
      ['FileButton','fileName','Name of export coordinate file:','Select export file',None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      #['SelectionListWithNone','useFixedAtomSerials','Use format atom serials:',allFormatsDict.keys(),None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],
     
   ['hExchProtection',
    'readHExchProtection',
    'writeHExchProtection',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','Proton exchange protection file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['hExchRate',
    'readHExchRate',
    'writeHExchRate',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','Proton exchange rates file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['shiftAnisotropy',
    'readShiftAnisotropy',
    'writeShiftAnisotropy',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','Chemical shift anisotropy file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['shifts',
    'readShifts',
    'writeShifts',
    {'auremol':  ['auremol.meta',1,0,         [],[]],
     'autoAssign':['autoassign.output',1,0,   [['CheckButton','getSeqFlag','Import sequence from file:',0,None]],
                                              []],
     'cns':      ['file.ppm',1,1,             [],[]],
     'cosmos':   ['cosmos.cod',0,1,           [],[]],
     'csi':      ['data.csi',1,1,             [],[]],
     'cyana':    ['',
                  'XEasy format',
                  'XEasy format',
                                            [],[]],
     'mars':     ['chemShifts.tab',1,1,     [],[]],
     'monte':    ['monte.cs',1,0,           [],[]],
     'nmrStar':  ['nmrStar.str',1,0,        [],[]],
     'nmrView':  ['ppm.out',1,1,            [],[]],
     'pipp':     ['pipp.shifts',1,1,        [['CheckButton','getSeqFlag','Import sequence from file:',0,None]],
                                            []],
     'pistachio':['molecule.csa',1,0,          [],[]],
     'pronto':   ['ss.report',1,0,          [],[]],
     'shiftx':   ['shifts.tab',1,1,         [],[]],
     'shifty':   ['shifts.tab',1,1,         [],[]],
     'sparky':   ['sparky.ppm',1,1,         [],[]],
     'talos':    ['project.tab',1,1,        [],[]],
     'targetedAcquisition':['TAassignment.ta', 1, 0, [], []],
     'xeasy':    ['xeasy.prot',1,1,         [],
                                            [['CheckButton','cyana21Naming','Use CYANA 2.1 atom names:',0,None]]]},
     [['FileButton','fileName','Chemical shift file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['SelectionListWithNone','measurementList','Add to existing shift list:',"project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList')",None],
      ['CheckButton','autoMergeShifts','Automatically merge shifts:',0,None],      
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [['SelectionList','measurementList','Select shift list to export:',"project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList')",None],
      ['FileButton','fileName','Name of export chemical shift file:','Select export file',None],
      ['MultiSelectionList','chains','Use only chains:','project.sortedMolSystems().sortedChains()',None],
      ['CheckButton','useOriginalResNames','Use original resonance names',0,None],
      ['CheckButton','compressResonances','Compress names in output',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['t1Relax',
    'readT1s',
    'writeT1s',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','T1 relaxation file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['t1RhoRelax',
    'readT1Rhos',
    'writeT1Rhos',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','T1 rho relaxation file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['t2Relax',
    'readT2s',
    'writeT2s',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','T2 relaxation file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['noe',
    'readNoes',
    'writeNoes',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','NOE file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['jCoupling',
    'readJCouplings',
    'writeJCouplings',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','J coupling file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['rdc',
    'readRdcs',
    'writeRdcs',
    {'nmrStar':  ['nmrStar.str',1,0,        [],[]]
    },
     [['FileButton','fileName','RDC file to import:','Select file',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['MultiSelectionList','experiments','Link to NMR experiments:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','measurementListName','List name:','',{'width': 20}],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [#See shifts for export info, if ever required... 
      ],
     None,
     None],

   ['peaks',
    'readPeaks',
    'writePeaks',
    {'ansig':    ['ansig.xpk',1,1,          [],
                                            [['MultiSelectionList','peakLists','Select peak lists to export:',"project.currentNmrProject.sortedExperiments().sortedDataSources().sortedPeakLists()",None],
                                             ['IntEntry','forceNumDim','Force number dimensions in output file (0 for automatic setting):',0,None]]],
     'auremol':  ['peaks.ml',1,0,           [],[]],
     'autoAssign':['autoAssign.pks',1,1,    [],[]],
     'ccpNmr':   ['ccpNmrAnalysis.tab',1,0, [],[]],
     'cyana':    ['spec-names.peaks',1,1,   [['CheckButton','forceHeight','Force volumes as height values:',0,None]],
                                            [['SelectionList','integrationMethod','Peak integration method:',['volume','height'],None]]],
     'felix':    ['peaks.txt',1,0,          [],[]],
     'marvin':   ['peaks.dat',1,0,          [['FloatEntry','probabilityCutoff','Set probability cutoff:',0.0,None],
                                             ['CheckButton','preservePeaks','Preserve existing peaks:',True,None],
                                             ['CheckButton','preserveAssignments','Preserve existing assignments:',True,None],
                                             #['CheckButton','removeDuplicateAssignments','Remove duplicate assignments:',True,None],
                                             ['StringEntry','defaultChainCode','Code of CCPN chain to link to:','A',None],
                                             ['StringEntry','defaultMolSystem','Name of CCPN molSystem to link to:','MS1',None],
                                             ['CheckButton','createNew','Create new peaks:',False,None],
                                             ['StringEntry','axisMapping','Order of marvin axis:','',None],
                                             ['CheckButton','useAppData','Use application data for peak mapping:',False,None],
                                             ['StringEntry','appDataKey','applicationName,keyword pair for searching appData:','',None]],
                                            []],
     'nmrDraw':  ['nmrDraw.tab',1,1,        [['CheckButton','usePointsInfo','Import width info as points:',0,None]],
                                            [['CheckButton','useAliasing','Export aliased chemical shifts:',1,None]]],
     'nmrView':  ['nmrView.xpk',1,1,        [['StringEntry','assignTagSep','Assignment separator:',' ',{'width': 1}],
                                             ['FileButton','processingParsFile','Spectrum parameter file:','Select file (optional)',{'component': 'processingPars'}]],
                                            [['CheckButton','writeMissingPeaks','Write missing peaks:',1,None],
                                             ['StringEntry','assignTagSep','Assignment separator:',' ',{'width': 1}]]],
     'pipp':     ['pipp.pck',1,0,           [],[]],
     'pronto':   ['pronto.xpk',1,0,         [],[]],
     'sparky':   ['sparky.peaks',1,1,       [],
                                            [['MultiSelectionList','volumeIntensityType','Intensity and/or volume export to file:',[volume_kw,height_kw],None]]],
     'xeasy':    ['xeasy.peaks',1,1,        [['CheckButton','cyanaFormat','Read as CYANA format:',0,None],
                                             ['CheckButton','forceHeight','Force volumes as height values:',0,None]],
                                            [['CheckButton','cyanaFormat','Write as CYANA format:',0,None],
                                             ['SelectionList','integrationMethod','Peak integration method:',['volume','height'],None]]]},
     [['FileButton','fileName','Peak list file to import:','Select file',None],
      ['SelectionListWithNone','dataSource','Use existing dataSource:','project.currentNmrProject.sortedExperiments().sortedDataSources()',None],
      ['SelectionListWithNone','peakList','Overwrite existing peaklist:','project.currentNmrProject.sortedExperiments().sortedDataSources().sortedPeakLists()',None],
      ['StringEntry','specName','Peak list name:','',{'width': 20}],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [['SelectionList','peakLists','Select peak list to export:',"project.currentNmrProject.sortedExperiments().sortedDataSources().sortedPeakLists()",None],
      ['FileButton','fileName','Name of export peak list file:','Select export file',None],
      ['CheckButton','forcePositiveVolumes','Force positive volumes/intensities in output:',0,None],
      ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
      ['CheckButton','compressResonances','Compress names in output:',1,None],
      ['CheckButton','writeAssignments','Include peak assignments in output:',1,None],
      ['CheckButton','useOriginalNumbers','Use original peak numbers:',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['peakAssignments',
    'readPeakAssignments',
    'writePeakAssignments',
    {'xeasy':    ['xeasy.assign',1,1,       [],[]]},
     [['FileButton','fileName','Peak assignment file to import:','Select file',None],
      ['SelectionListWithNone','peakList','Assignments for peaklist:','project.currentNmrProject.sortedExperiments().sortedDataSources().sortedPeakLists()',None],
      ['CheckButton','addMode','Add to existing assignments:',1,None],
      ['CheckButton','overwriteMode','Overwrite existing assignments:',0,None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [['SelectionList','peakList','Select peak list to export:',"project.currentNmrProject.sortedExperiments().sortedDataSources().sortedPeakLists()",None],
      ['FileButton','fileName','Name of export peak assignments file:','Select export file',None],
      ['CheckButton','useOriginalResNames','Use original resonance names',0,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],


   ['distanceConstraints',
    'readDistanceConstraints',
    'writeDistanceConstraints',
    {'amber':   ['amber.list',1,0,          [['FileButton','coordFile','Coordinate file (for atom names):','Select file (required)',None]],
                                            []],
     'aqua':    ['aqua.noe',1,1,            [],[]],
     'cns':     ['distances.tbl',1,1,       [],
                                            [['CheckButton','oneLinePerItem','Use one line per constraint item:',0,None],
                                             ['CheckButton','usePeakInfo','Use peak information:',1,None]]],
     'concoord':['distances.restr',1,1,     [],[]],
     'discover':['distances.restr',1,0,     [],[]],
     'dyana':   ['distances.upl',1,1,       [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None]],
                                            [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None]]],
     'cyana':   ['distances.upl',1,1,       [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None],
                                             ['SelectionList','version','Cyana version:',['2.0','2.1'],None]],
                                            [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None],
                                             ['SelectionList','version','Cyana version:',['2.0','2.1'],None]]],
     'mardigras':['distances.restr',1,0,    [],[]],
     'nmrStar': ['nmrStar.str',1,0,         [],[]],    },
     [['FileButton','fileName','Distance constraints file to import:','Select file',None],
      ['SelectionListWithNone','strucGen','Add to structure generation:','project.currentNmrProject.sortedStructureGenerations()',None],
      ['SelectionListWithNone','nmrConstraintStore','Add to NMR constraint store:','project.sortedNmrConstraintStores()',None],
      ['SelectionListWithNone','constraintList','Add to constraint list:',"project.currentNmrProject.sortedStructureGenerations().nmrConstraintStore.findAllConstraintLists(className = 'DistanceConstraintList')",None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [['SelectionList','constraintList','Select distance constraint list to export:',"project.currentNmrProject.sortedNmrConstraintStores().findAllConstraintLists(className = 'DistanceConstraintList')",None],
      ['FileButton','fileName','Name of export distance constraints file:','Select export file',None],
      ['CheckButton','sortForExport','Sort constraints for export:',0,None],
      ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
      ['CheckButton','compressResonances','Compress names in output:',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['hBondConstraints',
    'readHBondConstraints',
    'writeHBondConstraints',
    {'amber':   ['amber.list',1,0,          [['FileButton','coordFile','Coordinate file (for atom names):','Select file (required)',None]],
                                            []],
     'aqua':    ['aqua.noe',1,1,            [],[]],
     'cns':     ['hbonds.tbl',1,1,          [],[['CheckButton','oneLinePerItem','Use one line per constraint item:',0,None]]],
     'discover':['distances.restr',1,0,     [],[]],
     'dyana':   ['hbonds.upl',1,1,          [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None]],
                                            [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None]]],
     'cyana':   ['hbonds.upl',1,1,          [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None],
                                             ['SelectionList','version','Cyana version:',['2.0','2.1'],None]],
                                            [['FileButton','lowerLimitFile','Lower limit file:','Select file (optional)',None],
                                             ['SelectionList','version','Cyana version:',['2.0','2.1'],None]]],
     'nmrStar': ['nmrStar.str',1,0,         [],[]]},
     [['FileButton','fileName','H-bond constraints file to import:','Select file',None],
      ['SelectionListWithNone','strucGen','Add to structure generation:','project.currentNmrProject.sortedStructureGenerations()',None],
      ['SelectionListWithNone','nmrConstraintStore','Add to NMR constraint store:','project.sortedNmrConstraintStores()',None],
      ['SelectionListWithNone','constraintList','Add to constraint list:',"project.currentNmrProject.sortedStructureGenerations().nmrConstraintStore.findAllConstraintLists(className = 'HBondConstraintList')",None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],               
     [['SelectionList','constraintList','Select H-bond constraint list to export:',"project.currentNmrProject.sortedNmrConstraintStores().findAllConstraintLists(className = 'HBondConstraintList')",None],
      ['FileButton','fileName','Name of export H-bond constraints file:','Select export file',None],
      ['CheckButton','sortForExport','Sort constraints for export:',0,None],
      ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
      ['CheckButton','compressResonances','Compress names in output:',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['dihedralConstraints',
    'readDihedralConstraints',
    'writeDihedralConstraints',
    {'amber':   ['amber.list',1,0,          [['FileButton','coordFile','Coordinate file (for atom names):','Select file (required)',None]],
                                            []],
     'cns':     ['dihedrals.tbl',1,1,       [],[]],
     'discover':['dihedrals.restr',1,0,     [],[]],
     'dyana':   ['dihedrals.aco',1,1,       [],[]],
     'cyana':   ['dihedrals.aco',1,1,       [['SelectionList','version','Cyana version:',['2.0','2.1'],None]],
                                            [['SelectionList','version','Cyana version:',['2.0','2.1'],None]]],
     'nmrStar': ['nmrStar.str',1,0,         [],[]],
     'talos':   ['dihedral.tab',1,0,        [['MultiSelectionList','validClasses','Classes for valid constraint:',['Good','New','Warn','Dyn','Strong','Generous','None'],None],
                                             ['IntEntry','validCount','Count variable lower limit for valid constraint:',0,None],
                                             ['FloatEntry','multiplyDeviation','Multiplication factor to convert deviation in constraint range:',2.0,None]],
                                            []]},
     [['FileButton','fileName','Dihedral constraints file to import:','Select file',None],
      ['SelectionListWithNone','strucGen','Add to structure generation:','project.currentNmrProject.sortedStructureGenerations()',None],
      ['SelectionListWithNone','nmrConstraintStore','Add to NMR constraint store:','project.sortedNmrConstraintStores()',None],
      ['SelectionListWithNone','constraintList','Add to constraint list:',"project.currentNmrProject.sortedStructureGenerations().nmrConstraintStore.findAllConstraintLists(className = 'DihedralConstraintList')",None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],               
     [['SelectionList','constraintList','Select dihedral constraint list to export:',"project.currentNmrProject.sortedNmrConstraintStores().findAllConstraintLists(className = 'DihedralConstraintList')",None],
      ['FileButton','fileName','Name of export dihedral constraints file:','Select export file',None],
      ['CheckButton','sortForExport','Sort constraints for export:',0,None],
      ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
      ['CheckButton','compressResonances','Compress names in output:',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['rdcConstraints',
    'readRdcConstraints',
    'writeRdcConstraints',
    {'cns':     ['rdc.tbl',1,1,             [],
                                            [['IntEntry','refFrameResId','The residue code for the reference frame:',999,None],
                                             ['StringEntry','refFrameSegId','The segment Id for the reference frame:',defaultMolCode,{'width': 1}],
                                             ['CheckButton','writeWeight','Export value for restraint weight:',0,None],
                                             ['StringEntry','filterAtomNames','Only include atom names:','',{'width': 10}]]],
     'dyana':   ['rdcs',1,1,                [],
                                            [['CheckButton','singleAtom','Write out with single atom for each RDC:',0,None]]],
     'cyana':   ['rdcs',1,1,                [['SelectionList','version','Cyana version:',['2.0','2.1'],None]],
                                            [['SelectionList','version','Cyana version:',['2.0','2.1'],None],
                                             ['CheckButton','singleAtom','Write out with single atom for each RDC:',0,None]]],
     'module':  ['rdcs.mod',1,1,            [],[]],
     'pales':   ['rdcs.pal',1,1,            [],[]],
     'nmrStar': ['nmrStar.str',1,0,         [],[]]},
     [['FileButton','fileName','Rdc constraints file to import:','Select file',None],
      ['SelectionListWithNone','strucGen','Add to structure generation:','project.currentNmrProject.sortedStructureGenerations()',None],
      ['SelectionListWithNone','nmrConstraintStore','Add to NMR constraint store:','project.sortedNmrConstraintStores()',None],
      ['SelectionListWithNone','constraintList','Add to constraint list:',"project.currentNmrProject.sortedStructureGenerations().nmrConstraintStore.findAllConstraintLists(className = 'RdcConstraintList')",None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],               
     [['SelectionList','constraintList','Select RDC constraint list to export:',"project.currentNmrProject.sortedNmrConstraintStores().findAllConstraintLists(className = 'RdcConstraintList')",None],
      ['FileButton','fileName','Name of export RDC constraints file:','Select export file',None],
      ['CheckButton','negateValues','Negate exported coupling constant values:',0,None],
      ['CheckButton','sortForExport','Sort constraints for export:',0,None],
      ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
      ['CheckButton','compressResonances','Compress names in output:',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['csaConstraints',
    'readCsaConstraints',
    'writeCsaConstraints',
    {'cns':     ['csa.tbl',1,1,             [],
                                            [['IntEntry','refFrameResId','The residue code for the reference frame:',999,None],
                                             ['StringEntry','refFrameSegId','The segment Id for the reference frame:',defaultMolCode,{'width': 1}],
                                             ['CheckButton','writeWeight','Export value for restraint weight:',0,None]]],
                                             },
     [['FileButton','fileName','Chemical shift anisotropy (CSA) constraints file to import:','Select file',None],
      ['SelectionListWithNone','strucGen','Add to structure generation:','project.currentNmrProject.sortedStructureGenerations()',None],
      ['SelectionListWithNone','nmrConstraintStore','Add to NMR constraint store:','project.sortedNmrConstraintStores()',None],
      ['SelectionListWithNone','constraintList','Add to constraint list:',"project.currentNmrProject.sortedStructureGenerations().nmrConstraintStore.findAllConstraintLists(className = 'ChemShiftConstraintList')",None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],               
     [['SelectionList','constraintList','Select chemical shift anisotropy (CSA) constraint list to export:',"project.currentNmrProject.sortedNmrConstraintStores().findAllConstraintLists(className = 'ChemShiftConstraintList')",None],
      ['FileButton','fileName','Name of export CSA constraints file:','Select export file',None],
      ['CheckButton','sortForExport','Sort constraints for export:',0,None],
      ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
      ['CheckButton','compressResonances','Compress names in output:',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['jCouplingConstraints',
    'readJCouplingConstraints',
    'writeJCouplingConstraints',
    {'cns':     ['jcouplings.tbl',1,0,       [],[]]},
     [['FileButton','fileName','J coupling constraints file to import:','Select file',None],
      ['SelectionListWithNone','strucGen','Add to structure generation:','project.currentNmrProject.sortedStructureGenerations()',None],
      ['SelectionListWithNone','nmrConstraintStore','Add to NMR constraint store:','project.sortedNmrConstraintStores()',None],
      ['SelectionListWithNone','constraintList','Add to constraint list:',"project.currentNmrProject.sortedStructureGenerations().nmrConstraintStore.findAllConstraintLists(className = 'JCouplingConstraintList')",None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],               
     [['SelectionList','constraintList','Select J coupling constraint list to export:',"project.currentNmrProject.sortedNmrConstraintStores().findAllConstraintLists(className = 'JCouplingConstraintList')",None],
      ['FileButton','fileName','Name of export J coupling constraints file:','Select export file',None],
      ['CheckButton','sortForExport','Sort constraints for export:',0,None],
      ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
      ['CheckButton','compressResonances','Compress names in output:',1,None],
      ['CheckButton','useOriginalData','Use original data:',0,None],
      ['CheckButton','resetMapping','Reset chain mapping:',0,None],
      ['CheckButton','useCcpnChainInfo','Use CCPN chain and sequence codes:',0,None],
      ['SelectionListWithNone','forceNamingSystemName','Use naming system:',namingSystemNames,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],

   ['acqPars',
    'readAcqPars',
    'writeAcqPars',
    {'bruker':  ['acqus',1,0,               [['CheckButton','switchRefDimInfo','Switch acqu2s/3s info:',1,None]],
                                            []],
     'varian':  ['procpar',1,0,             [],[]]},
     [['FileButton','fileName','Acquisition parameters file to import:','Select file',None],
      ['SelectionListWithNone','experiment','Use experiment:','project.currentNmrProject.sortedExperiments()',None],
      ['SelectionListWithNone','fidDataSource','Use FID dataSource:',"project.currentNmrProject.sortedExperiments().findAllDataSources(dataType = 'FID')",None],
      #['SelectionListWithNone','freqDataSource','Use FREQ dataSource:',"project.currentNmrProject.experiments.findAllDataSources(dataType = 'processed')",None],
      ['CheckButton','overwritePars','Overwrite existing parameters:',0,None],               
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],               
     [],
     None,
     ["AcqProcParsEditPopup(self.guiParent, experiment = formatObject.experiment, newFidDataSources = ('nmrPipe','azara'))"]],

   ['processingPars',
    'readProcPars',
    'writeProcPars',
    {'azara':  ['',1,1,                     [],[]],
     'bruker':   ['',1,0,                   [],[]],
     'nmrPipe':  ['',0,1,                   [],[]],
     'nmrView':  ['',1,0,                   [],[]]},
     [['FileButton','fileName','Processing parameters file to import:','Select file',None],
      ['SelectionListWithNone','experiment','Use experiment:','project.currentNmrProject.sortedExperiments()',None],
      ['StringEntry','expName','Name for experiment:','',None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],               
     [['SelectionList','freqDataSource','Select FREQ data source to export:',"project.currentNmrProject.sortedExperiments().findAllDataSources(dataType = 'processed')",None],
      ['FileButton','fileName','Name of export processing parameters file:','Select export file',None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     ["AcqProcParsEditPopup(self.guiParent, experiment = addKeywds['freqDataSource'].experiment, newFidDataSources = ('nmrPipe','azara'))",
      'return'],
     None],

   ['processingScript',
    'readProcScript',
    'writeProcScript',
    {'azara':  ['',0,1,                     [],[]],
     'nmrPipe':  ['',0,1,                   [],[]]},
     [],
     [['SelectionList','freqDataSource','Select FREQ data source to export:',"project.currentNmrProject.sortedExperiments().findAllDataSources(dataType = 'processed')",None],
      ['FileButton','fileName','Name of export processing script file:','Select export file',None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     ["AcqProcParsEditPopup(self.guiParent, experiment = addKeywds['freqDataSource'].experiment, newFidDataSources = ('nmrPipe','azara'))",
      'return'],
     None],

   ['peopleAndCitations',
    'readPeopleAndCitations',
    '',
    {'pdb':     ['code.pdb',1,0,                [],[]]},
     [['FileButton','fileName','File with people/citation information to import:','Select file',None],
      ['CheckButton','useExistingPersons','Match to Persons already in project:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [],
     None,
     None],

   ['assignment',
    'readAssignment',
    '',
    {'mars':     ['assignment_AA.out',1,0,                [['SelectionList','reliabilityCutoff','Reliability cutoff:',['H','M','L'],None]],
                                                          []]},
     [['FileButton','fileName','File with assignments to import:','Select file',None],
      ['SelectionList','chain','Select chain to link imported data to:','project.sortedMolSystems().sortedChains()',None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [],
     None,
     None],
   
   ['project',
    'readProject',
    'writeProject',
    {'ansig':      ['ansig.ctr',1,0,           [['FileButton','peakFile','Select associated peak file:',"Select peak file",{'component': 'peaks'}]],
                                                []],
     'autoAssign': ['table.aat',1,1,           [],
                                               [['SelectionList','chain','Select chain to analyze:',"project.sortedMolSystems().sortedChains()",None]]],
     'cosmos':     ['cosmos.cod',0,1,          [],
                                               [['SelectionList','projectType','Project file (COO for coord/shifts, COD for restraints/shifts):',['COO','COD'],None],
                                                ['MultiSelectionList','models','Select models to export:',"project.sortedMolSystems().sortedStructureEnsembles().sortedModels()",None],
                                                ['SelectionList','shiftList','Select shift list to export:',"project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList')",None],
                                                ['SelectionList','constraintList','Select distance constraint list to export:',"project.currentNmrProject.sortedNmrConstraintStores().findAllConstraintLists(className = 'DistanceConstraintList')",None]]],
     'mars':       ['mars.inp',0,1,            [],
                                               [['SelectionList','chain','Select chain to analyze:',"project.sortedMolSystems().sortedChains()",None],
                                                ['DirButton','inputDirectory','Directory to write the input files to:','Select input directory',None],
                                                ['SelectionList','shiftList','Select shift list to use in the MARS input:',"project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList')",None],
                                                ['CheckButton','useFixedConnection','Write out fixed connections file:',0,None],
                                                ['CheckButton','useFixedAssignments','Write out fixed assignments file:',0,None]
                                                ]],
     'nmrStar':    ['bmrb.str',1,1,            [],
                                               [['SelectionList','entry','Select entry to export:',"project.currentNmrEntryStore.sortedEntries()",None],
                                                ['CheckButton','useOriginalData','Use original data:',0,None]]],
     'nmrView':    ['nmrView.str',1,0,         [],[]], 
     'sparky':     ['sparky.proj',1,1,         [['CheckButton','ignoreShifts','Ignore reading of chemical shifts:',0,None]],
                                               [['SelectionList','entry','Select entry to export:',"project.currentNmrEntryStore.sortedEntries()",None],
                                                ['SelectionList','molSystem','Select molecular system to export:',"project.sortedMolSystems()",None],
                                                ['CheckButton','useOriginalResNames','Use original resonance names:',0,None],
                                                ['CheckButton','useOriginalNumbers','Use original peak numbers:',1,None]]],
     'talos':      ['project.tab',1,1,         [],
                                               [['SelectionList','chain','Select chain to export:',"project.sortedMolSystems().sortedChains()",None],
                                                ['SelectionList','shiftList','Select shift list to export:',"project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList')",None],
                                                ['CheckButton','useOriginalData','Use original data:',0,None]]]},
     [['FileButton','fileName','Project file to import:','Select file',None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [['FileButton','fileName','Name of export project file:','Select export file',None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     None,
     None],
     
   #
   # Combined reading/writing
   #

   ['shiftsAndPeaks',
    'readShiftsAndPeaks',
    '',
    {'nmrView':  ['',1,0,            [],[]],
     'pronto':   ['',1,0,            [],[]],
     'sparky':   ['',1,0,            [],[]],
     'xeasy':    ['',1,0,            [['CheckButton','cyanaFormat','Read as CYANA format:',0,None],
                                      ['FileButton','peakAssignmentFile','Peak assignment file to import:','Select file',{'component': 'peakAssignments'}]],
                                     []]},
     [['FileButton','shiftFile','Chemical shift file to import:','Select file',{'component': 'shifts'}],
      ['FileButton','peakFile','Peak file to import:','Select file',{'component': 'peaks'}],
      ['SelectionListWithNone','dataSource','Use existing dataSource:','project.currentNmrProject.sortedExperiments().sortedDataSources()',None],
      ['CheckButton','resonanceMatching','Match existing resonances:',1,None],
      ['CheckButton','minimalPrompts','Minimal interaction:',0,None]],
     [],
     None,
     None],

]

##############
# File types #
##############

fileTypeDict = {
  
  'sequence':   {'generic':    ['*.seq*','*.pdb','*.m*l*'],
                 'ansig':      ['*.seq'],
                 'auremol':    ['*.seq'],
                 'autoAssign': ['*.aat'],
                 'charmm':     ['*.crd'],
                 'cns':        ['*.pdb'],
                 'csi':        ['data.*'],
                 'dyana':      ['*.seq'],
                 'cyana':      ['*.seq'],
                 'fasta':      ['*.seq'],
                 'gromos':     ['*.gro'],
                 'mars':       ['*.seq'],
                 'mol2':       ['*.mol2','*.ml2'],
                 'molmol':     ['*.pdb'],
                 'monte':      ['*.seq'],
                 'nmrStar':    ['*.str'],
                 'nmrView':    ['*.seq'],
                 'pdb':        ['*.pdb'],
                 'pistachio':  ['*.csa'],
                 'pronto':     ['*.report'],
                 'pseudoPdb':  ['*.pdb'],
                 'shiftx':     ['*.tab'],
                 'sparky':     ['*.seq'],
                 'talos':      ['*.tab'],
                 'xeasy':      ['*.seq']},

  'chemComps':  {'generic': ['*.pdb','*.m*l*'],
                 'auremol': ['*.comp'],
                 'mol':     ['*.mol'],
                 'mol2':    ['*.mol2','*.ml2'],
                 'pdb':     ['*.pdb']},

  'coordinates':{'generic':   ['*.pdb'],
                 'auremol':   ['*.pdb'],
                 'charmm':    ['*.crd'],
                 'cns':       ['*.pdb'],
                 'cosmos':    ['*.coo'],
                 'cyana':     ['*.cor'],
                 'dyana':     ['*.cor'],
                 'gromos':    ['*.gro'],
                 'molmol':    ['*.pdb'],
                 'pdb':       ['*.pdb'],
                 'pseudoPdb': ['*.pdb']},
  
  'hExchProtection': {'generic': ['*.str'],
                      'nmrStar': ['*.str']},

  'hExchRates':  {'generic': ['*.str'],
                  'nmrStar': ['*.str']},

  'shiftAnisotropy': {'generic': ['*.str'],
                      'nmrStar': ['*.str']},

  'shifts':         {'generic': ['*.ppm','*.prot'],
                     'auremol': ['*.meta'],
                     'autoAssign': ['*.out'],
                     'cns':     ['*.ppm'],
                     'cosmos':    ['*.cod'],
                     'csi':     ['data.*'],
                     'mars':    ['*.tab'],
                     'monte':   ['*.cs'],
                     'nmrStar': ['*.str'],
                     'nmrView': ['*.out'],
                     'pipp':    ['*.shifts'],
                     'pistachio':  ['*.csa'],
                     'pronto':  ['*.report'],
                     'shiftx':  ['*.tab'],
                     'shifty':  ['*.tab'],
                     'sparky':  ['*.ppm'],
                     'talos':   ['*.tab'],
                     'targetedAcquisition' : ['*.ta'],
                     'xeasy':   ['*.prot']},

  't1Relax':        {'generic': ['*.str'],
                     'nmrStar': ['*.str']},

  't1RhoRelax':     {'generic': ['*.str'],
                     'nmrStar': ['*.str']},

  't2Relax':        {'generic': ['*.str'],
                     'nmrStar': ['*.str']},

  'noe':            {'generic': ['*.str'],
                     'nmrStar': ['*.str']},

  'peaks':      {'generic': ['*.xpk','*.peaks','*.tbl','*.tab'],
                 'ansig':   ['*.xpk','*.export'],
                 'auremol': ['*.ml'],
                 'autoAssign': ['*.pks'],
                 'cyana':   ['*.peaks'],
                 'felix':   ['*.txt'],
                 'nmrDraw': ['*.tab','*.tbl'],
                 'nmrView': ['*.xpk'],
                 'pipp':    ['*.pck','*.PCK'],
                 'pronto':  ['*.xpk'],
                 'sparky':  ['*.peaks'],
                 'xeasy':   ['*.peaks']},

  'peakAssignments':  {
  
                 'xeasy': ['*.assign']},

  'distanceConstraints':   {
  
                 'generic': ['*.tbl','*.restr','*.upl','*.lol'],
                 'amber':   ['*.list'],
                 'aqua':    ['*.noe'],
                 'cns':     ['*.tbl'],
                 'concoord':['*.restr'],
                 'discover':['*.restr'],
                 'dyana':   ['*.upl','*.lol'],
                 'cyana':   ['*.upl','*.lol'],
                 'mardigras':['*.restr','*.rstrnt'],
                 'nmrStar': ['*.str']},

  'hBondConstraints':   {
  
                 'generic': ['*.tbl','*.upl','*.lol'],
                 'amber':   ['*.list'],
                 'aqua':    ['*.noe'],
                 'cns':     ['*.tbl'],
                 'discover':['*.restr'],
                 'dyana':   ['*.upl','*.lol'],
                 'cyana':   ['*.upl','*.lol'],
                 'nmrStar': ['*.str']},

  'dihedralConstraints':   {
  
                 'generic': ['*.tbl','*dihed*'],
                 'amber':   ['*.list'],
                 'cns':     ['*.tbl'],
                 'discover':['*.restr'],
                 'dyana':   ['*.aco'],
                 'cyana':   ['*.aco'],
                 'nmrStar': ['*.str'],
                 'talos':   ['*.tab']},

  'rdcConstraints':   {
  
                 'generic': ['*.tbl','*rdc*'],
                 'cns':     ['*.tbl'],
                 'dyana':   ['*rdc*'],
                 'cyana':   ['*rdc*'],
                 'pales':   ['*.pal','*.txt'],
                 'module':  ['*.mod'],
                 'nmrStar': ['*.str']},

  'csaConstraints':   {
  
                 'generic': ['*.tbl'],
                 'cns':     ['*.tbl']},

  'jCouplingConstraints':   {
  
                 'generic': ['*.tbl'],
                 'cns':     ['*.tbl']},

  'acqPars':    {'generic': ['acqu*','procpar'],
                 'bruker':  ['acqu*'],
                 'varian':  ['procpar']},

  'acqPars':    {'generic': ['acqu*','procpar'],
                 'bruker':  ['acqu*'],
                 'varian':  ['procpar']},

  'processingPars': {

                 'azara':   ['*.par'],
                 'bruker':  ['proc*'],
                 'nmrPipe': ['*.com'],
                 'nmrView': ['*.par']},

  'processingScript': {

                 'azara':   ['*.scr'],
                 'nmrPipe': ['*.com']},
  
  'peopleAndCitations':  {'pdb':     ['*.pdb']},
  
  'assignment':   {'mars':     ['*.out']},

  'project':    {'ansig':      ['*.ctr'],
                 'autoAssign': ['*.aat'],
                 'cosmos':     ['*.coo'],
                 'mars':       ['*.inp'],
                 'nmrStar':    ['*.str'],
                 'nmrView':    ['*.str'],
                 'sparky':     ['*.proj','*.save'],
                 'talos':      ['*.tab']},

  'mapping':    {ccpNmr_kw: ['*.map']},

  'csdep':      {'bmrb': ['*.csdep']},

  'FID':        {'varian': ['fid'],'bruker': ['ser','fid']},

  #'conv':       {'nmrPipe': ['*.com']},

  'processing': {'nmrPipe': ['*.com'], 'azara': ['*.scr','*.par']},
  
  #'par': {'azara': ['*.par']}
  
  'text':       {'generic': ['*.txt','*.text']}

}

##########################################################################
# Conversion from old experiment names to new reference experiment setup #
##########################################################################

oldToNewExpDict = {

  #
  # 2 dimensional HSQC type experiments
  # -----------------------------------

  'hsqc_HC.hc': ('H[C]','H[C]'),
  'hsqc_HN.hn': ('H[N]','H[N]'),

  # Tocsy spectra
  # -------------

  'tocsy.hh': ('H_H.relayed','H_H.relayed'),

  'tocsy_HcCh.hc': ('HC_CH.relayed','hC_cH.relayed'),
  'tocsy_HCCH.hhc': ('HC_CH.relayed','HC_cH.relayed'),

  'tocsy_CcconH.hc': ('HC_CCONH.relayed','hC_cconH.relayed'),
  'tocsy_CccoNH.hcn': ('HC_CCONH.relayed','hC_ccoNH.relayed'),
  'tocsy_HccoNH.hhn': ('HC_CCONH.relayed','Hc_ccoNH.relayed'),

  'tocsy_CcNH.hcn': ('HC_CNH.relayed','hC_cNH.relayed'),
  'tocsy_HcNH.hhn': ('HC_CNH.relayed','Hc_cNH.relayed'),

  # Cosy spectra
  # ------------

  'cosy.hh': ('H_H.COSY','H_H.COSY'),

  # Backbone (including beta) assignment spectra. Grouped by complement
  # SP is info about same or previous residue, P is previous residue only
  # ---------------------------------------------------------------------

  'HNCO.hcn': ('H[N[CO]]','H[N[CO]]'),
  'HnCO.hc': ('H[N[CO]]','H[n[CO]]'),

  'HNcoCA.hcn': ('H[N[CO[CA]]]','H[N[co[CA]]]'),
  'HncoCA.hc': ('H[N[CO[CA]]]','H[n[co[CA]]]'),

  'HNcaCO.hcn': ('H[N[CA[CO]]]','H[N[CA[CO]]] CA splits CO'),
  'HncaCO.hc': ('H[N[CA[CO]]]','H_H.COSY'),

  'HNCA.hcn': ('H[N[CA]]','H[N[CA]]'),
  'HnCA.hc': ('H[N[CA]]','H[n[CA]]'),

  'haCANH.hcn': ('HNCAHA','haCANH'),
  'HAcaNH.hhn': ('HNCAHA','HAcaNH'),

  'haCAcoNH.hcn': ('HCACONH','hCAcoNH'), 
  'HAcacoNH.hhn': ('HCACONH','HcacoNH'), 
  
  'hbCBCAcoNH.hcn': ('H{CA|Cca}CONH','h{CA|Cca}coNH'), 

  'hbCBCANH.hcn': ('H{CA|Cca}NH','h{CA|Cca}NH'),
  
  'HNcoCACB.hcn': ('H[N[CO[{CA|ca[C]}]]]','H[N[co[{CA|ca[C]}]]]'),    
  'HncoCACB.hc': ('H[N[CO[{CA|ca[C]}]]]','H[n[co[{CA|ca[C]}]]]'),
   # TODO are two below correct?
  'HNcocaCB.hcn': ('H[N[CO[{CA|ca[C]}]]]','H[N[co[{CA|ca[C]}]]]'), 
  'HncocaCB.hc': ('H[N[CO[{CA|ca[C]}]]]','H[n[co[{CA|ca[C]}]]]'),

  'HNcaCB.hcn': ('H[N[CA[CB]]]','H[N[ca[CB]]]'),
  'HncaCB.hc': ('H[N[CA[CB]]]','H[n[ca[CB]]]'),
  
  'HNCACB.hcn': ('H[N[{CA|ca[Cali]}]]','H[N[{CA|ca[Cali]}]]'),
  'HnCACB.hc': ('H[N[{CA|ca[Cali]}]]','H[n[{CA|ca[Cali]}]]'),  
  
  #trosy_HnCO.hc   # SP. Same as above, but TROSY version.
  #trosy_HncaCO.hc # P. Same as above, but TROSY version.
  #trosy_HncoCA.hc # P. Same as above, but TROSY version.
  #trosy_HnCA.hc   # SP. Same as above, but TROSY version.
  #trosy_HnCACB.hc # SP. Same as above, but TROSY version.
  #trosy_HncoCACB.hc # P. Same as above, but TROSY version.
  #trosy_HncaCB.hc # SP. Same as above, but TROSY version.
  #trosy_HncocaCB.hc # P. Same as above, but TROSY version.
  
  #
  # Noesy spectra
  #
  
  'noesy.hh': ('H_H.through-space','H_H.through-space'),
  'noesy_hsqc_HCH.hhc': ('H[C]_H.through-space','H_H[C].through-space'),
  'noesy_hsqc_HNH.hhn': ('H[N]_H.through-space','H_H[N].through-space'),
  'noesy_hsqc_HXH.hhx': ('H[{N|C}]_H.through-space','H_H[{N|C}].through-space'),
  'noesy_hsqc_HCCH.hhcc': ('H[C]_H[C].through-space','H[C]_H[C].through-space'),
  'noesy_hsqc_HNNH.hhnn': ('H[N]_H[N].through-space','H[N]_H[N].through-space'),
  'noesy_hsqc_HCNH.hhcn': ('H[N]_H[C].through-space','H[C]_H[N].through-space'),
  'noesy_hsqc_HNCH.hhnc': ('H[N]_H[C].through-space','H[N]_H[C].through-space')
  
}

########################
# Software definitions #
########################

softwareDefinitions = {

     # This one for Analysis!
     'ccpNmr':{'version': 'any',
      
                  'tasks': ['Spectrum display','Spectrum analysis'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Volume',
                          'Intensity'],

                     'Peak picking':
                 
                         ['Automatic and manual']
 
                     }
                  },

      # TODO: if making changes here, also do for cyana!
      'xeasy':  {'version': 'any',
      
                 'tasks': ['Spectrum display','Spectrum analysis'],
                 
                 'methods': {
                 
                     'Peak integration':

                         ['Denk integration',
                          'Rectangular integration',
                          'Elliptical integration',
                          'Maximum integration',
                          'Automatic integration'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                    }

                  },

      'cyana':  {'version': 'any',
      
                 'tasks': ['Spectrum display','Spectrum analysis'],
                 
                 'methods': {
                 
                     'Peak integration':

                         ['Denk integration',
                          'Rectangular integration',
                          'Elliptical integration',
                          'Maximum integration',
                          'Automatic integration'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                    }

                  },
                
       'nmrView':{'version': 'any',
      
                  'tasks': ['Spectrum display','Spectrum analysis'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Volume',
                          'Intensity'],

                     'Peak picking':
                 
                         ['Automatic and manual']
 
                     }
                  },
  
       'auremol':{'version': 'any',
      
                  'tasks': ['Spectrum display','Spectrum analysis','Structure calculation'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Volume',
                          'Intensity'],

                     'Peak picking':
                 
                         ['Automatic and manual']
 
                     }
                  },

       'azara':  {'version': 'any',
      
                  'tasks': ['Spectrum display','Spectrum analysis'],
                 
                  'methods': {}
                 
                 },
                 
       'ansig':  {'version': 'any',
      
                  'tasks': ['Spectrum display','Spectrum analysis'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Volume'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                     }
                 },
                  
       'nmrPipe':{'version': 'any',
      
                  'tasks': ['Spectrum processing'],
                 
                  'methods': {}
                  
                },

       'pipp':  {'version': 'any',
      
                 'tasks': ['Spectrum analysis'],
                 
                 'methods': {
                 
                     'Peak integration':

                         ['Intensity'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                     }
                 },

       'pronto': {'version': 'any',
      
                  'tasks': ['Spectrum analysis'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Intensity'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                     }
                 },
 
       'unknown': {'version': 'any',
      
                  'tasks': ['Spectrum analysis'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Intensity'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                     }
                 },

       'sparky': {'version': 'any',
      
                  'tasks': ['Spectrum display','Spectrum analysis'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Volume (Lorentzian)',
                          'Volume (Gaussian)',
                          'Volume (Sum over box)',
                          'Volume (Sum over ellipse)',
                          'Volume (Unknown method)',
                          'Intensity (data)',
                          'Intensity (fitted)'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                                        }
                 },
                              
       'nmrDraw':{'version': 'any',
      
                  'tasks': ['Spectrum display','Spectrum analysis'],
                 
                  'methods': {
                 
                     'Peak integration':

                         ['Volume',
                          'Intensity'],

                     'Peak picking':
                 
                         ['Automatic and manual']
                     
                                        }
                 },
                 
       'felix':{'version': 'any',
      
                'tasks': ['Spectrum processing','Spectrum display','Spectrum analysis'],
                 
                'methods': {
                
                   'Peak integration':

                       ['Volume',
                        'Intensity'],

                   'Peak picking':
                 
                       ['Automatic and manual']
 
                                      }
               }
  }


######################################
# Processing/acquisition definitions #
######################################

varianAcqOrdInfo = {3: [['d3,d2,ph2,ph (normal)','d3,d2,ph,ph2'],
                        [('dim3','dim2','phase2','phase1'),('dim3','dim2','phase1','phase2')]],
                    4: [['d4,d3,d2,ph3,ph2,ph (normal)','d4,d3,d2,ph,ph2,ph3'],
                        [('dim4','dim3','dim2','phase3','phase2','phase1'),('dim4','dim3','dim2','phase1','phase2','phase3')]]
                   }

"""
   
xeasySpinSysAtoms = ['N','HN','CA','HA','HA2','CB','HB2','HB3','QB','CG',
                     'HG2','HG3','QG','CD','HD2','HD3','QD','CE','HE2',
                     'HE3','QE','H?1','H?2','H?3','H?4','H?5','HA1','C?1',
                     'C?2','C?3','C?4','N?1','C']   

"""
