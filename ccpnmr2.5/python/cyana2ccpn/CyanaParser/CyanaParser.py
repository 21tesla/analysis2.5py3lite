#!/usr/bin/env python
#
#
"""
CyanaParser setup

run as:
   python CyanaParser.py


Data structure
    CyanaParser:
            ->sequence:list of SequenceRecord instances
            ->atoms:list of AtomRecord instances
            ->molecule:MoleculeRecord instance
            ->convention:cyana convention string (CYANA or CYANA2)
            ->resonances:list of Resonance instances
            ->peakLists:list of PeakList instances
            ->seqFile:xeasy seq file
            ->protFile:xeasy prot file
            ->pdbFile:cyana pdb file
            ->modelCount

    MoleculeRecord:
        ->chains:list of ChainRecord instances

    ChainRecord:
        ->chainId
        ->molecule:MoleculeRecord instance
        ->residues:list of ResidueRecord instances

    ResidueRecord:
        ->chainId
        ->sequenceId
        ->residueId
        ->typeIdentifier:ResidueDef instance (if defined)
        ->chain:ChainRecord instance
        ->atoms:list of AtomRecord instances

    AtomRecord:
        ->chainId
        ->sequenceId
        ->residueId
        ->atomId
        ->resonance:Resonance instance
        ->xeasyIndex:original Xeasy index from prot file (if known)
        ->pdbRecords:list of PDBrecord instances for atom (as obtained from PyMMlib)
        ->residue:ResidueRecord instance

    Resonance:
        ->value
        ->error
        ->stereo
        ->atom:AtomRecord instance

    PeakList:list of Peak instances

    Peak:
        ->dimension
        ->positions:list of floats
        ->height:NTvalue instance
        ->resonances:list of Resonance instances
        ->xeasyIndex:original Xeasy index from prot file (if known)

    PDBRecord:
        ->atom:AtomRecord instance
        ->model
        ->x
        ->y
        ->z
        ->occupancy
        ->tempFactor


Adapted from PluginCode xeasy.py and cyana.py files

"""

import sys
#sys.path.append('.')

# from cing.core.constants import * #@UnusedWildImport
import cing.Libs.NTutils as ntu
from cing.Libs.AwkLike import AwkLike
#from cing.Libs.disk import path
from cing.core.classes import Peak
from cing.core.classes import PeakList
from cing.core.classes import DistanceRestraintList
from cing.core.classes import DistanceRestraint
from cing.core.classes import DihedralRestraintList
from cing.core.classes import DihedralRestraint
from cing.core.molecule import Resonance
from cing.core.database import translateAtomName
from cing.Libs import PyMMLib
from cing.Libs.pdb import moveFirstDigitToEnd
from cing.Libs.pdb import MatchGame
from cyana2ccpn.classes4 import RDCRestraint
from cyana2ccpn.classes4 import RDCRestraintList
from cyana2ccpn.classes4 import ChemicalShiftRestraint
from cyana2ccpn.classes4 import ChemicalShiftRestraintList
import json

CYANA_CHAIN_ID = 'A'
CYANA_NON_RESIDUES = ['PL','LL2','link']
CYANA2     = 'CYANA2'
CYANA     = 'CYANA'
X_AXIS = 0
INTERNAL_0 = 'INTERNAL_0'   # INTERNAL_0 is the first convention used: was based upon DYANA/CYANA1.x convention (Gly has HA1/2)
INTERNAL_1 = 'INTERNAL_1'   # INTERNAL_1 is the second convention used: IUPAC for IUPAC defined atoms, CYANA2 for non-IUPAC atoms
INTERNAL   = INTERNAL_1


class MoleculeRecord( ntu.NTtree ):
    """MoleculeRecord class; root of the chain,residue,atoms records tree
    """
    def __init__(self):
        ntu.NTtree.__init__(self, 'molecule' )
    #end def

    @property
    def chains(self):
        return self._children

    def toTuple(self):
        return( ('MoleculeRecord',self.name) )
    #end def

    def __str__(self):
        return ntu.sprintf('<MoleculeRecord %s>', self.name)
    #end def

    def __repr__(self):
        return ntu.sprintf('MoleculeRecord()')
#end class

class ChainRecord( ntu.NTtree ):
    """ChainRecord class; 'defined' by tuple (chainId)
    """
    def __init__(self, chainId):
        ntu.NTtree.__init__(self, str(chainId) )
        self.chainId = chainId
    #end def

    @property
    def residues(self):
        return self._children

    @property
    def molecule(self):
        return self._parent

    def toTuple(self):
        return( ('ChainRecord',self.chainId) )
    #end def

    def __str__(self):
        return ntu.sprintf('<ChainRecord %s>', self.chainId)
    #end def

    def __repr__(self):
        return ntu.sprintf('ChainRecord(%s)', self.chainId)
#end class


class ResidueRecord( ntu.NTtree ):
    """ResidueRecord class; 'defined' by tuple (chainId, sequenceId, residueId)
    """
    def __init__(self, chainId, sequenceId, residueId):
        ntu.NTtree.__init__(self, str(residueId)+str(sequenceId) )
        self.chainId = chainId
        self.sequenceId = sequenceId
        self.residueId = residueId
        self.typeIdentifier = None
    #end def

    @property
    def atoms(self):
        return self._children

    @property
    def chain(self):
        return self._parent

    def toTuple(self):
        return( ('ResidueRecord',self.chainId,self.sequenceId,self.residueId) )
    #end def

    def __str__(self):
        return ntu.sprintf('<ResidueRecord %s.%s%s>', self.chainId, self.residueId, self.sequenceId)
    #end def

    def __repr__(self):
        return ntu.sprintf('ResidueRecord(%s,%s,%s)', self.chainId, self.sequenceId, self.residueId)
    #end def
#end class

class AtomRecord( ntu.NTtree ):
    """AtomRecord class; 'defined' by tuple (chainId, sequenceId, residueId, atomId)
    """
    nextId=0
    def __init__(self, chainId, sequenceId, residueId, atomId):
        ntu.NTtree.__init__(self, str(atomId))
        self.chainId    = chainId
        self.sequenceId = sequenceId
        self.residueId  = residueId
        self.atomId     = atomId
        self.id         = AtomRecord.nextId # DihedralRestraint needs an unique id attribute
        AtomRecord.nextId          += 1
        self.typeIdentifier = None
        self.resonance  = None   #link to resonance object
        self.xeasyIndex = None   #original Xeasy index from prot file (if known)
        self.pdbRecords = ntu.NTlist() # list of PDBrecord instances as obtained from PyMMlib
    #end def

    @property
    def residue(self):
        return self._parent

    def toTuple(self):
        return( ('AtomRecord',self.chainId,self.sequenceId,self.residueId,self.atomId) )
    #end def

    def __str__(self):
        return ntu.sprintf('<AtomRecord %s.%s%s.%s>', self.chainId, self.residueId, self.sequenceId, self.atomId)
    #end def

    def __repr__(self):
        return ntu.sprintf('AtomRecord(%s,%s,%s,%s)', self.chainId, self.sequenceId, self.residueId, self.atomId)
    #end def
#end class


def findTupleInDict( myDict, *args):
    """Convenience methods to find a tuple defined by args.
    returns item if present or None otherwise.
    test for empty dicts.
    """
    if myDict == None or len(myDict) == 0:
        return None
    if len(args) == 0:
        return None
    t = tuple(args)
    if myDict.has_key(t):
        return myDict[t]
    else:
        return None
#end def


class CyanaParser(dict):
    """
    Wrapper around the Cing Xeasy/Cyana import routines
    """

    def __init__(self):
        dict.__init__(self)                   # dict contains all sequence and atom tuples as keys and points to its instances;
                                              # i.e. can be used as a lookup for validity
        self.sequence       = ntu.NTlist()    # list of SequenceRecord instances
        self.atoms          = ntu.NTlist()    # list of AtomRecord instances
        self.convention     = None            # CYANA or CYANA2 (3) convention
        self.resonances     = ntu.NTlist()    # list of Resonance instances
        self.peakLists      = ntu.NTlist()    # list of PeakList instances
        self.distanceRestraintLists = ntu.NTlist()    # list of DistanceRestraintList instances
        self.dihedralRestraintLists = ntu.NTlist()    # list of DistanceRestraintList instances
        self.chemicalShiftRestraints = ntu.NTlist()
        self.violationLists = ntu.NTlist()
        self.rdcRestraintLists = ntu.NTlist()    # list of RDCRestraintList instances
        self.pdbRecords     = None            # PDB records as returned by PYMMLib
        self.modelCount     = 0               # number od models in the pdb file
        self.molecule       = None            # MoleculeRecord instance

        self.rootPath       = None            # rootPath to directory with cyana files
        self.seqFile        = None
        self.protFile       = None
        self.pdbFile        = None
        self.peakFiles      = None
        self.cyaFiles       = None
        self.uplFiles       = None
        self.txtFiles       = None
        self.lolFiles       = None
        self.acoFiles       = None
        self.rdcFiles       = None

        self._seqDict       = None            # interim dict for seqFile, protFile, peaksFile parsing; sequenceId mapping
        self._protDict      = None            # interim dict for seqFile, protFile, peaksFile parsing; AtomRecord.xeasyIndex mapping
        self._atomDict      = None            # interim dict for stereo file or .upl / .lol file parsing;
                                              #         (AtomRecord.sequenceId,AtomRecord) mapping
    #end def

    def parse(self,
              cyanaPath = None,
              seqFile   = None,
              finalProtFile  = None,
              originalProtFile  = None,
              pdbFile   = None,
              cyaFiles  = ntu.NTlist(),
              peakFiles = ntu.NTlist(),
              xpkFiles  = ntu.NTlist(),
              noaFiles  = ntu.NTlist(),
              txtFiles  = ntu.NTlist(),
              uplFiles  = ntu.NTlist(),
              lolFiles  = ntu.NTlist(),
              acoFiles  = ntu.NTlist(),
              rdcFiles  = ntu.NTlist(),
             ):
        """Extract data from either cyanaPath directory, relying on file extensions to decide on content
        or using explicit file definitions
        Only allows for one .seq file, one .prot file and one .pdb file
        return self on None on error
        """
        self.seqFile    = seqFile
        self.finalProtFile   = finalProtFile
        self.originalProtFile   = originalProtFile
        self.pdbFile    = pdbFile
        self.cyaFiles   = cyaFiles
        self.peakFiles  = peakFiles
        self.xpkFiles   = xpkFiles
        self.noaFiles   = noaFiles
        self.txtFiles   = txtFiles
        self.uplFiles   = uplFiles
        self.lolFiles   = lolFiles
        self.acoFiles   = acoFiles
        self.rdcFiles   = rdcFiles

        if cyanaPath != None:
            # extract files from cyanaPath using extentions
            self.rootPath = path(cyanaPath)
            if not self.rootPath.exists() or not self.rootPath.isdir():
                ntu.nTerror('CyanaParser.parse: invalid directory %s', cyanaPath)
                return None
            #end if
            ntu.nTmessage('Converting files from directory %s', cyanaPath)

            # find files from extension
            for f in self.rootPath:
                dir,base,ext = f.split3()
                if ext == '.seq':
                    self.seqFile = f
                elif ext == '.prot':
                    self.protFile = f
                elif ext == '.pdb':
                    self.pdbFile = f
                elif ext == '.cya':
                    self.cyaFiles.append(f)
                elif ext == '.peaks':
                    self.peakFiles.append(f)
                elif ext == '.xpk':
                    self.xpkFiles.append(f)
                elif ext == '.upl':
                    self.uplFiles.append(f)
                elif ext == '.noa':
                    self.noaFiles.append(f)
                elif ext == '.txt':
                    self.txtFiles.append(f)
                elif ext == '.lol':
                    self.lolFiles.append(f)
                elif ext == '.aco':
                    self.acoFiles.append(f)
                #~ elif ext == '.rdc':
                    #~ self.rdcFiles.append(f)
                #end if
            #end for
        else:
            self.rootPath = 'Undefined'
        #end if


        self._atomDict = {}
        # parse seqFile and finalProtFile to assemble the sequence, atoms and resonances lists
        if self.seqFile != None and self.finalProtFile != None:
            self.parseSeqFile( self.seqFile, init=True )
            self.parseProtFiles( self.originalProtFile, self.finalProtFile, init=True )
        #end if



        if self.pdbFile:
            self.parsePdbFile( self.pdbFile )
        #end if

        #set convention and try to match
        self.convention = self.autoDetectCyanaConvention()
        self._matchResiduesAndAtoms()

        # parse all other files
        for f in self.cyaFiles:
            self.parseCyanaStereoFile(f)
        #end for

        for f in self.peakFiles:
            peakList = self.parsePeakFile(f)
            if peakList != None:
                self.peakLists.append(peakList)
            #end if
        #end for
        for f in self.xpkFiles:
            peakList = self.parseXpkFile(f)
            if peakList != None:
                self.peakLists.append(peakList)
            #end if
        #end for
        ccpnConfigFile = open("Properties.ccpn.json")
        ccpnConfig = json.load(ccpnConfigFile)
        violationLimit = float(ccpnConfig["violationLimit"])
        for f in self.noaFiles:
                  violationList = self.parseNoaFile(f, violationLimit)
                  distanceList = self.getDistanceRestraintsFromNoa(f)
                  if violationList != None:
                      self.violationLists.append(violationList)
                  if distanceList != None:
                    self.distanceRestraintLists.append(distanceList)
                  #end if
              #end for

        for f in self.uplFiles:
            restraintList = self.parseUplFile(f)
            if restraintList != None:
                self.distanceRestraintLists.append(restraintList)
            #end if
        #end for

        for f in self.lolFiles:
            pass
        #end for

        for f in self.acoFiles:
            restraintList = self.parseAcoFile(f)
            if restraintList != None:
                self.dihedralRestraintLists.append(restraintList)
            #end if
        #end for

        for f in self.rdcFiles:
            restraintList = self.parseRdcFile(f)
            if restraintList != None:
                self.rdcRestraintLists.append(restraintList)

        #end for

        return self
    #end def

    def findTuple(self, *args):
        """find item define by args as tuple in self
        return item or None if not present
        """
        return findTupleInDict(self,*args)
    #end def

    def appendObject(self, object):
        """
        Append ChainRecord, ResidueRecord or AtomRecord to self
        Recursively add chain and molecule if needed
        Maintain linkage and referencing in self
        return object or None on error
        """
        if object == None:
            return None

        if self.has_key(object.toTuple()):
            ntu.nTerror('CyanaParser.appendObject: object %s already present', object)
            return None
        #end if

        if isinstance(object, MoleculeRecord):
            self.molecule = object
            self[object.toTuple()] = object

        elif isinstance(object, ChainRecord):
            if self.molecule == None:
                self.appendObject( MoleculeRecord() )
            #end if
            self.molecule.addChild2( object )
            self[object.toTuple()] = object

        elif isinstance(object, ResidueRecord):
            chain = self.findTuple('ChainRecord', object.chainId)
            if chain == None:
                chain = self.appendObject( ChainRecord(object.chainId) )
            chain.addChild2(object)
            self.sequence.append(object)
            self[object.toTuple()] = object

        elif isinstance(object, AtomRecord):
            self.atoms.append(object)
            self[object.toTuple()] = object
            self._atomDict[(object.sequenceId,object.atomId)] = object # alternative references that are handy/needed.
            res = self.findTuple('ResidueRecord', object.chainId, object.sequenceId, object.residueId)
            if res == None:
                ntu.nTerror('CyanaParser.appendObject: unable to identify residue for %s', object)
            else:
                res.addChild2( object )
            #end if
        else:
            ntu.nTerror('CyanaParser.appendObject: invalid opject type %s', type(object))
            return None
        #end if
        return object
    #end def

    def parseSeqFile( self, seqFile, init=True):

        """parse the seqFile, fill _seqDict and if init==True also sequence
        return self or None on error
        """
        self._seqDict = {}
        if init: self.sequence = ntu.NTlist()
        sequenceId = 1
        resCount = 0
        error = False
        for line in AwkLike( seqFile, commentString='#' ):
            if (not line.isEmpty() and not line.isComment( '#')):
                residueId = line[1]
                if ( residueId in CYANA_NON_RESIDUES ):        # skip the bloody CYANA non-residue stuff
                    pass

                else:
                    if (line.NF > 1):
                        sequenceId = line.int(2)
                        if sequenceId == None:
                            error = True
                            ntu.nTerror('CyanaParser.parseSeqFile: invalid sequenceId "%s" on line %d (%s)',
                                        line[2], line.NR, line[0] )
                        #end if
                    #endif
                    self._seqDict[ sequenceId ] = residueId # store original 'convention' name
                    if init:
                        self.appendObject( ResidueRecord(CYANA_CHAIN_ID,sequenceId, residueId) )
                    sequenceId += 1
                    resCount += 1
                #end if
            #end if
        #end for
        self.seqFile = seqFile
        ntu.nTmessage('==> Parsed %d residues from %s', resCount, self.seqFile)

        if error:
            return None
        else:
            return self
    #end def

    def parseProtFiles( self, originalProtFile, finalProtFile, init=True):

        """parse the protFile, fill _protDict and if init==True also atoms
        requires _seqDict to be initialised (using parseSeqFile method)
        return self or None on error
        """
        if self._seqDict == None or len(self._seqDict) == 0:
            ntu.nTerror('CyanaParser.parseProtFile: uninitialised seqDict; call parseSeqFile first')
            return None
        #end if

        self._protDict = {}
        if init: self.atoms = ntu.NTlist()
        self.protCount = 0
        error = False

        for line in  AwkLike( finalProtFile, commentString='#' ):
            if line.NF == 5:
                # Xeasy/Cyana atom index
                index     = line.int( 1 )
                shift     = line.float( 2 )
                error     = line.float( 3 )
                atomId    = line[4]
                sequenceId = line.int( 5 )
                if sequenceId not in self._seqDict:
                    ntu.nTwarning( 'CyanaParser.parseFinalProtFile: undefined sequenceId %d in "%s:%d" (%s)' % (
                                    sequenceId, protFile, f.NR, f[0]))
                    error = True
                else:
                    residueId = self._seqDict[sequenceId]

                    if init:
                        atm = AtomRecord(CYANA_CHAIN_ID,sequenceId,residueId,atomId)
                        atm.xeasyIndex = index
                        self.appendObject( atm )
                    else:
                        atm = self.findTuple('AtomRecord', CYANA_CHAIN_ID,sequenceId,residueId,atomId)
                    #end if

                    if atm == None:
                        ntu.nTwarning( 'CyanaParser.parseFinalProtFile: cannot define atom in "%s:%d" (%s)',
                                        protFile, f.NR, f[0]
                                     )
                        error = True
                        continue
                    #end if
                    if shift != 999.000:
                        reson = Resonance( atm, shift, error )
                        reson.stereo = False
                        self.resonances.append(reson )
                    else:
                        reson = None
                    #end if
                    atm.resonance = reson # reverse linkage
                    self._protDict[ index ] = atm
                    self.protCount += 1
                #end if
            #end if
        #end for



    	for line in  AwkLike( originalProtFile, commentString='#' ):
            if line.NF == 5:
                # Xeasy/Cyana atom index
                index     = line.int( 1 )
                shift     = line.float( 2 )
                error     = line.float( 3 )
                atomId    = line[4]
                sequenceId = line.int( 5 )
                residueId = self._seqDict[sequenceId]

                if init:
                        atom = AtomRecord(CYANA_CHAIN_ID,sequenceId,residueId,atomId)
                        atom.xeasyIndex = index
                        #self.appendObject( atm )
                else:
                        atom = self.findTuple('AtomRecord', CYANA_CHAIN_ID,sequenceId,residueId,atomId)
                    #end if

                if atom == None:
                        ntu.nTwarning( 'CyanaParser.parseFinalProtFile: cannot define atom in "%s:%d" (%s)',
                                        originalProtFile, f.NR, f[0]
                                     )
                        error = True
                        continue
                if shift != 999.000:
                			reson = Resonance( atom, shift, error )
                			if reson in self.resonances:

                				atom.resonance = reson
                				chemicalShift = ChemicalShiftRestraint( atom=atom, value=shift, error=error )
                				self.chemicalShiftRestraints.append( chemicalShift )


        self.finalProtFile = finalProtFile
        self.originalProtFile= originalProtFile
        ntu.nTmessage('==> Parsed %d atoms from %s', self.protCount, self.protFile)
        if error:
            return None
        else:
            return self
    #end def

    def _parseXpkLine( self, line):
        """Parse line taking whitespace and first {} as field delimeters.
        return list with fields or None on error
        """
        import string

        result = []

        l = len(line)
        i = 0
        tokenStart = False
        parsingToken = False
        curlyToken = False
        while i<l: # loop over all characters; has to be dynamic

            isWhiteSpace = line[i] in string.whitespace

            #print i, l, parsingToken, curlyToken, isWhiteSpace, c[i]

            if not parsingToken:
                if isWhiteSpace:
                    i += 1
                elif line[i] == '{':
                    parsingToken = True
                    curlyToken = True
                    i += 1
                    tokenStart = i
                elif line[i] == '}':
                    i += 1
                else:
                    parsingToken = True
                    curlyToken = False
                    tokenStart = i
                    i += 1
            else:
                if curlyToken:
                    if line[i] == '}':
                        result.append(line[tokenStart:i])
                        parsingToken = False
                        curlyToken = False
                        i += 1
                    elif line[i] == '{':
                        i += 1
                        tokenStart = i
                    else:
                        i += 1
                    #end if
                else:
                    if isWhiteSpace:
                        result.append(line[tokenStart:i])
                        parsingToken = False
                        curlyToken = False
                        i += 1
                    else:
                        i += 1
                #end if
            #end if
        #end while
        print len(result),'>>', result
        return result
    #end def

    def parseXpkFile( self, xpkFile)   :
        """Read nmrView xpk peak file
           returns a PeakList instance or None on error

label dataset sw sf
H h C
None
{10000.00 } {5147.74 } {4001.28 }
{500.1300 } {500.1300 } {125.7570 }
 H.L  H.P  H.W  H.B  H.E  H.J  H.U  h.L  h.P  h.W  h.B  h.E  h.J  h.U  C.L  C.P  C.W  C.B  C.E  C.J  C.U  vol  int  stat  comment  flag0
1  {?}   7.101   0.050   0.050   ?   0.000   {?}   {32.HB2}   2.361   0.050   0.050   ?   0.000   {?}   {?}   130.550   0.050   0.050   ?   0.000   {?}  25782257.12500 1011583.00000 1 {?} 0
2  {?}   5.736   0.050   0.050   ?   0.000   {?}   {?}   6.954   0.050   0.050   ?   0.000   {?}   {?}   51.989   0.050   0.050   ?   0.000   {?}  21689452.09375 1544707.87500 1 {?} 0
3  {13.HB2}   1.982   0.050   0.050   ?   0.000   {?}   {?}   0.989   0.050   0.050   ?   0.000   {?}   {?}   41.360   0.050   0.050   ?   0.000   {?}  17728599.21484 870235.00000 1 {?} 0
4  {?}   3.820   0.050   0.050   ?   0.000   {?}   {?}   8.082   0.050   0.050   ?   0.000   {?}   {?}   63.644   0.050   0.050   ?   0.000   {?}  48441156.75000 4535237.50000 1 {?} 0
5  {45.HA}   4.434   0.050   0.050   ?   0.000   {?}   {?}   7.465   0.050   0.050   ?   0.000   {?}   {?}   57.981   0.050   0.050   ?   0.000   {?}  75860595.21875 5855013.00000 1 {?} 0
6  {45.HA}   4.434   0.050   0.050   ?   0.000   {?}   {46.H}   7.619   0.050   0.050   ?   0.000   {?}   {?}   57.996   0.050   0.050   ?   0.000   {?}  47893949.96875 3641987.00000 1 {?} 0
7  {?}   4.118   0.050   0.050   ?   0.000   {?}   {?}   2.971   0.050   0.050   ?   0.000   {?}   {?}   55.169   0.050   0.050   ?   0.000   {?}  14386389.29688 1090827.12500 1 {?} 0
8  {?}   3.816   0.050   0.050   ?   0.000   {?}   {{82.H 83.H 85.H 88.H}}   8.520   0.050   0.050   ?   0.000   {?}   {?}   63.102   0.050   0.050   ?   0.000   {?}  43355361.82812 2821826.00000 1 {?} 0

problem resides in the multiple assignments that alters format
        """

        error = False

        _path,name,_ext = path( xpkFile ).split3()
        peaks = PeakList( name=name, status='keep' )

        dimension = 0
        fieldNames = []
        header = {}
        # f stands for field.
        for line in  AwkLike( xpkFile ):

            #print line, line.NF, dimension, fieldNames, header

            # header info in the top
            if (line.NR == 1 ):                  # first line
                fieldNames = line[1:]
            elif (line.NR <= len(fieldNames)+1): # header lines
                #print line.NR,fieldNames, header
                header[fieldNames[line.NR-2]] = line[1:]
            elif (line.NR == len(fieldNames)+2): # field names line
                header[ 'fields' ] = line[1:]
                print 'header>', header

            # 'regular' lines with data
            elif (not line.isComment('#') and line.NF > 0):

                # try to determine dimension from header info
                if header.has_key('label'):
                    dimension = len( header['label'] )
                #end if

                if not dimension:
                    ntu.nTerror('CyanaParser.parseXpkFile: invalid dimensionality in file "%s" (line %d, "%s")'%(
                                 xpkFile, line.NR, line[0]))
                    return None
                #end if

                # we need to reparse the line because of syntax
                tokens = self._parseXpkLine(line[0])
                # and update 'line'
                for i,token in enumerate(tokens):
                    line[i+1] = token
                #end for
                line.NF = len(tokens)

                # pick up relevant elements
                cur = 1

                # preserve the peak id
                peakId = line.int( cur )
                if (peakId == None):
                    return None
                cur += 2


                peakpos = []
                for _i in range(X_AXIS, dimension):
                    p = line.float( cur )
                    if (p == None):
                        return None
                    peakpos.append( p )
                    cur += 4
                #end if
                cur -= 2 # skip two fields
                heightField = line[cur].split('*')
                height = float(heightField[-1])
                if height == None:
                    return None
                #cur += 1
                #heightError = line.float( cur )
                #if heightError == None:
                #    return None
                #cur += 1
                #
                resonances = []
                error = False
                cur = 2 # skip two fields
                for _i in range(X_AXIS, dimension):
                    print _i
#                     aIndex = line.int( cur )
#                     if aIndex == None:
                #        return None
                #    cur += 1
                #    # index 0 means unassigned according to Xeasy convention
                #    if aIndex == 0:
                #        resonances.append( None )
                #    else:
                #        if not aIndex in self._protDict:
                #            ntu.nTerror('CyanaParser.parsePeakFile: invalid atom id %d on line %d (%s)',
                #                         aIndex, line.NR, line[0]
                #                       )
                #            error = True
                #            break
                #        else:
                #            resonances.append(self._protDict[ aIndex].resonance)
                #        #end if
                #    #end if
                ##end for

                #if not error:
                #    peak = Peak( dimension=dimension,
                #                 positions=peakpos,
                #                 height=height, heightError=heightError,
                #                 resonances = resonances,
                #                )
                #    # store original peak id
                #    peak.xeasyIndex = peakId
                #    peaks.append( peak )
                ##end if
            #end if
        #end for

        peaks.xpkFile = xpkFile
        ntu.nTmessage('==> Parsed %d peaks from %s', len(peaks), xpkFile )
        #end if

        return peaks
    #end def

    def parsePeakFile( self, peakFile)   :
        """Read Xeasy peak file
           returns a PeakList instance or None on error

           JFD: description of XEASY peak list format:
  43   1.760   3.143 1 T          0.000e+00  0.00e+00 -   0 2260 2587 0
  46   1.649   4.432 1 T          1.035e+05  0.00e+00 r   0 2583 2257 0
   ^ peak id                      ^ height
       ^ chemical shifts                     ^ height dev   ^ resonance ids
                     ^ ?                              ^ ?             ^ ?
                       ^ ?                                ^ ?

        resonance id is zero for unassigned.

        30.01.14 SPS: amended to also parse peak files exported by CYANA, which often contain
        multiple assignments per peak
        """

        name,_ext = peakFile.split('.peaks')
        peaks = PeakList( name=name, status='keep' )

        dimension = 0
        # f stands for field.
        for line in  AwkLike( peakFile ):
            if (line.NR == 1 and line.NF == 5):
                dimension = line.int(5)

            elif (not line.isComment('#') ):

                if not dimension:
                    ntu.nTerror('CyanaParser.parsePeakFile: invalid dimensionality in file "%s" (line %d, "%s")'%(
                                 peakFile, line.NR, line[0]))
                    return None
                #end if

                resonances = []
                ambiguous = False

                if line.NF <= dimension+6:
                	# Ambiguous assignment for this peak
									ambiguous = True
# 									print line
									for i in range(dimension):
											aIndex = line.int(i+1)
											if aIndex == 0:
													resonances.append( None )
											else:
												if not aIndex in self._protDict:
															ntu.nTerror('CyanaParser.parsePeakFile: invalid atom id %d on line %d (%s)',
																					 aIndex, line.NR, line[0]
																				 )
															error = True
															break
												else:
															resonances.append(self._protDict[ aIndex].resonance)
												#end if
											#end if
									#end for

                else:
									cur = 1

									# preserve the Xeasy peak id
									peakId = line.int( cur )
									if (peakId == None):
											return None
									cur += 1

									peakpos = []
									for _i in range(X_AXIS, dimension):
											p = line.float( cur )
											if (p == None):
													return None
											peakpos.append( p )
											cur += 1
									#end for

									cur += 2 # skip two fields
									height = line.float( cur )
									if height == None:
											return None
									cur += 1
									heightError = line.float( cur )
									if heightError == None:
											return None
									cur += 1

									resonances = []
									error = False
									cur += 2 # skip two fields
									for _i in range(X_AXIS, dimension):
											aIndex = line.int( cur )
											if aIndex == None:
													return None
											cur += 1
											# index 0 means unassigned according to Xeasy convention
											if aIndex == 0:
													resonances.append( None )
											else:
													if not aIndex in self._protDict:
															ntu.nTerror('CyanaParser.parsePeakFile: invalid atom id %d on line %d (%s)',
																					 aIndex, line.NR, line[0]
																				 )
															error = True

															break
													else:
															resonances.append(self._protDict[ aIndex].resonance)
													#end if
											#end if
									#end for
								#endif

                if not error:
                	if not ambiguous:
                	  peak = Peak( dimension=dimension,
                                 positions=peakpos,
                                 height=height, heightError=heightError,
                                 resonances = resonances,
                                )
                	  # store original peak id
                	  peak.xeasyIndex = peakId
                	  peaks.append( peak )
                	else:
                	  peak = peaks()
                	  for r in resonances:
                  		peak.resonances.append(r)
                #end if
            #end if
        #end for

        peaks.peakFile = peakFile
        ntu.nTmessage('==> Parsed %d peaks from %s', len(peaks), peakFile )
        #end if

        return peaks
    #end def

    def autoDetectCyanaConvention(self):
        'Returns None on error or CYANA or CYANA2'
        if len(self.atoms) <= 0:
            return None

        countMap = ntu.CountMap()
        countMap['HN'] = 0 # CYANA 1
        countMap['H'] = 0
        for a in self.atoms:
            atmName = a.atomId
            if not (atmName == 'HN' or atmName == 'H'):
                continue
            #end if
            countMap.increaseCount(atmName, 1)
        #end for
        overallCount = countMap.overallCount()
        if not overallCount:
            ntu.nTmessage("Assuming default CYANA2 convention because no amide protons counted at all.")
            return CYANA2
        #end if
        #msg = "Amide protons counted of CYANA 1/2 (HN/H) are: %s/%s." % ( countMap['HN'], countMap['H'])
        if countMap['HN'] > countMap['H']:
            ntu.nTmessage("==> Autodetected CYANA1.x convention")
            return CYANA
        else:
            ntu.nTmessage("==> Autodetected CYANA2.x / CYANA3.x convention")
            return CYANA2
    #end def

    def _matchResiduesAndAtoms(self):
        """Match residues to Cing NTdb entries
        """
        match = MatchGame(convention=self.convention, patchAtomNames = False, skipWaters = False, allowNonStandardResidue = False)
        errorCount = 0
        for res in self.sequence:
            # tmp add resName
            res.resName = res.residueId
            db=match.matchResidue2Cing(res)
            del res.resName
            if db != None:
                res.typeIdentifier = db
            else:
                ntu.nTerror('CyanaParser._matchResiduesAndAtoms: unable to determine type of residue %s', res)
                errorCount += 1
        #end for
        for atm in self.atoms:
            db=match.matchAtom2Cing(atm)
            if db != None:
                atm.typeIdentifier = db
            else:
                ntu.nTerror('CyanaParser._matchResiduesAndAtoms: unable to determine type of atom %s', res)
                errorCount += 1
        #end for
        if errorCount == 0:
            ntu.nTmessage("==> Successfully matched all residues and atoms")
        else:
            ntu.nTmessage("==> Matched residues and atoms but encountered %d errors", errorCount)
    #end def

    def _setStereo(self, sequenceId, atomId):
        """set stereo assignment of specific atom
        return AtomRecord or None on error
        """

        atm = findTupleInDict(self._atomDict, sequenceId, atomId)
        if atm == None:
            return None
        elif atm.resonance == None:
            return None
        else:
            atm.resonance.stereo = True
        return atm
    #end def

    def parseCyanaStereoFile( self, stereoFileName ):
        """Parse stereo assignments from CYANA stereo.cya type file
           return self or None on error.

    CYANA stereo file:

    var info echo
    echo:=off
    info:=none
    atom stereo "HB2  HB3   :509"   # GLU-
    atom stereo "QG1  QG2   :511"   # VAL
    atom stereo "HB2  HB3   :513"   # HIS
    atom stereo "QG1  QG2   :514"   # VAL
    atom stereo "HG2  HG3   :516"   # GLU-
    atom stereo "HA1  HA2   :519"   # GLY

        """
        count = 0
        for line in AwkLike( stereoFileName, minNF=5 ):
            if line[1] == 'atom' and line[2] == 'stereo':
                sequenceId = int (line[5].strip('"').strip(":") )
                for i in [3,4]:
                    atomId = line[i].strip('"')
                    atm = None
                    atm = self._setStereo( sequenceId, atomId)
                    if atm == None:
                        ntu.nTerror('CyanaParser.parseCyanaStereoFile: setting atom %s on line %d (%s)', atomId, line.NR, line[0] )
                        break
                    #end if
                    count += 1

                    # Val, Leu methyls: Carbon implicit in CYANA defs
                    implicit = [('VAL','QG1','CG1'),
                                ('VAL','QG2','CG2'),
                                ('LEU','QD1','CD1'),
                                ('LEU','QD2','CD2'),
                               ]
                    for resType, atmH, atmC in implicit:
                        if resType == atm.residueId and atmH == atm.atomId:
                            self._setStereo( atm.sequenceId, atmC )
                            count += 1
                        #end if
                    #end for
                #end for
            #end if
        #end for
        ntu.nTmessage('==> Derived %d stereo assignments from %s', count, stereoFileName )
        return self
    #end def


    def getDistanceRestraintsFromNoa(self, noaFile):

      maxErrorCount = 50
      errorCount = 0
      name,_ext = noaFile.split('.noa')
      if name == 'cycle7':
        listName = 'final'
      else:
        listName = name

      result = DistanceRestraintList( name=listName, status='keep')

      for line in AwkLike ( noaFile ):
        if line.NF >=7 and line[1] == 'Peak' and line[line.index('ppm;')+1] != 'diagonal):':
            peakId  = line[2]
            peakListName = line[4].split('.peaks')[0]
            upper = line[line.index('ppm;')+1]
            assignments = int(line.next()[4])
            for i in range(assignments):
              nextLine = line.next()
              if 'Violated' not in nextLine and '+' in nextLine:
                atmIdxList = [[3,1],[7,5]]
                atmList = []
                if '!' in nextLine:
                  i = nextLine.index('!')
                  del nextLine[i]
                if '*' in nextLine:
                  i = nextLine.index('*')
                  del nextLine[i]
                for atmIdx in atmIdxList:
                  sequenceId = int(nextLine[atmIdx[0]])
                  atomId = nextLine[atmIdx[1]]
                  atm = findTupleInDict(self._atomDict, sequenceId, atomId)
                  if not atm:
                      if errorCount <= maxErrorCount:
                          ntu.nTerror('CyanaParser.parseUplFile: Failed to decode for sequenceId %d atomId %s; line: %s', sequenceId, atomId, line[0] )
                      if errorCount == maxErrorCount+1:
                          ntu.nTerror("And so on")
                      errorCount += 1
                      continue
                  #end if
                  atmList.append( atm )
              #end for
                if len(atmList) != 2:
                  continue
                atm1 = atmList[0]
                atm2 = atmList[1]
                r = DistanceRestraint( atomPairs= [(atm1,atm2)], lower=0.0, upper=upper, peakList=peakListName+'-'+name, peak = peakId)
                result.append(r)
      return result

    def parseNoaFile (self, noaFile, violationLimit):
      """
      Peak 1740 from 1308PAR_6ms_SSB2p5.peaks (21.79, 21.14 ppm; 5.50 A):
      0 out of 4 assignments used, quality = 0.00:
        CG2   THR    7 - CG1   VAL    5  far    0    25   0   -  8.7-9.9
        CG2   THR    7 - CG1   VAL   70  far    0    22   0   -  9.1-12.1
        CG2   THR    7 - CG2   VAL   70  far    0    23   0   -  9.3-12.0
        CG2   THR   14 - CG1   VAL    5  far    0    25   0   -  9.6-10.7
      Violated in 20 structures by 2.66 A.


      """

      maxErrorCount = 50
      errorCount = 0
      name,_ext = noaFile.split('.noa')
      result = DistanceRestraintList( name=name, status='keep')
      violatedPeaks = {}
      for line in AwkLike ( noaFile ):
         if line.NF >=7 and line[1] == 'Peak':
           violatedPeaks[line[4].split('.peaks')[0]] = []
      for line in AwkLike ( noaFile ):
         if line.NF >=7 and line[1] == 'Peak':
            peakId  = line[2]
            peakListName = line[4].split('.peaks')[0]
         if line.NF > 5 and line[1] == 'Violated' and line[3] != "0" and float(line[6]) >= float(violationLimit):
            violatedPeaks[peakListName].append(peakId)
      for line in AwkLike ( noaFile ):
        if line.NF >=7 and line[1] == 'Peak' and line[line.index('ppm;')+1] != 'diagonal):':
            peakId  = line[2]
            peakListName = line[4].split('.peaks')[0]
            upper = line[line.index('ppm;')+1]
            # if len(line) == 12:
            #   upper = line[10]
            # elif len(line) == 11:
            #   upper = line[9]
            # elif len(line) == 10:
            #   upper = line[8]
            # print len(line), upper
            if peakId in violatedPeaks[peakListName]:
              assignments = int(line.next()[4])
              for i in range(assignments):
                nextLine = line.next()
                if 'Violated' not in nextLine:
                  atmIdxList = [[3,1],[7,5]]
                  atmList = []
                  if '!' in nextLine:
                    i = nextLine.index('!')
                    del nextLine[i]
                  if '*' in nextLine:
                    i = nextLine.index('*')
                    del nextLine[i]
                  for atmIdx in atmIdxList:
                    sequenceId = int(nextLine[atmIdx[0]])
                    atomId = nextLine[atmIdx[1]]
                    atm = findTupleInDict(self._atomDict, sequenceId, atomId)
                    if not atm:
                        if errorCount <= maxErrorCount:
                            ntu.nTerror('CyanaParser.parseUplFile: Failed to decode for sequenceId %d atomId %s; line: %s', sequenceId, atomId, line[0] )
                        if errorCount == maxErrorCount+1:
                            ntu.nTerror("And so on")
                        errorCount += 1
                        continue
                    #end if
                    atmList.append( atm )
                #end for
                  if len(atmList) != 2:
                    continue
                  atm1 = atmList[0]
                  atm2 = atmList[1]

                  r = DistanceRestraint( atomPairs= [(atm1,atm2)], lower=0.0, upper=upper, peakList=peakListName+'-'+name, peak = peakId)
                  result.append(r)
      return result


    def parseUplFile( self, uplFile, lower = 0.0 ):
        """
        Read Cyana upl file
        return a DistanceRestraintList or None on error
        """
        maxErrorCount = 50
        errorCount = 0

        if self._atomDict == None or len(self._atomDict) == 0:
            return None

        name,_ext = uplFile.split('.upl')
        result = DistanceRestraintList( name=name, status='keep')
        noQfRestraints = 0
        for line in AwkLike( uplFile, minNF=7 ):

          if '#QF' in line:

            # find the atoms
            atmIdxList = [[1,3],[4,6]]
            atmList = []
            for atmIdx in atmIdxList:

                sequenceId = line.int(atmIdx[0])
                atomId = line[atmIdx[1]]
                atm = findTupleInDict(self._atomDict, sequenceId, atomId)
                if not atm:
                    if errorCount <= maxErrorCount:
                        ntu.nTerror('CyanaParser.parseUplFile: Failed to decode for sequenceId %d atomId %s; line: %s', sequenceId, atomId, line[0] )
                    if errorCount == maxErrorCount+1:
                        ntu.nTerror("And so on")
                    errorCount += 1
                    continue
                #end if
                atmList.append( atm )
            #end for
            if len(atmList) != 2:
                continue
            # Unpack convenience variables.
            atm1 = atmList[0]
            atm2 = atmList[1]

            upper = line.float(7)
            # ambiguous restraint, should be append to last one
            if upper == 0:
                result().appendPair( (atm1,atm2) )
                continue
            if not upper:
                ntu.nTerror("CyanaParser.parseUplFile: Skipping line without valid upper bound on line: [" + line[0]+']')
                continue

            r = DistanceRestraint( atomPairs= [(atm1,atm2)], lower=lower, upper=upper )
            result.append( r )
            # also store the Candid info if present
            # print line
            if line.NF >= 9:
                r.peak = line.int( 9 )
                for peakList in self.peakLists:
                  for peak in peakList:
                    if peak.xeasyIndex == line.int(9):
                      peakAtoms = []
                      if None not in peak.resonances:
                        for resonance in peak.resonances:
                          peakAtoms.append(resonance.atom)
                      if atm1 in peakAtoms and atm2 in peakAtoms:
                        r.peakList=peakList.name
            if line.NF >= 11:
                r.SUP = line.float( 11 )
            if line.NF >= 13:
                r.QF = line.float( 13 )
          else:
            noQfRestraints+=1



        #end for
        if errorCount:
            ntu.nTerror("CyanaParser.parseUplFile: Found number of errors importing upl file: %s" % errorCount)

        result.uplFile = uplFile

        ntu.nTmessage('==> Imported upl file: %s from %s', result, uplFile )
        ntu.nTmessage('==> Number of Low Quality Restraints skipped: %s', noQfRestraints)
        return result
    #end def
    def parseRdcFile( self, rdcFile, lower = 0.0 ):
        """
        Read Cyana upl file
        return an RDCRestraintList or None on error
        """
        maxErrorCount = 50
        errorCount = 0

        if self._atomDict == None or len(self._atomDict) == 0:
            return None

        _dir,name,_ext = path( rdcFile ).split3()
        result = RDCRestraintList( name=name, status='keep')

        for line in AwkLike( rdcFile, commentString="#", minNF=7 ):

            # find the atoms
            atmIdxList = [[1,3],[4,6]]
            atmList = []
            for atmIdx in atmIdxList:
                #~
                sequenceId = line.int(atmIdx[0])
                atomId = line[atmIdx[1]]
                atm = findTupleInDict(self._atomDict, sequenceId, atomId)
                if not atm:
                    if errorCount <= maxErrorCount:
                        ntu.nTerror('CyanaParser.parseRdcFile: Failed to decode for sequenceId %d atomId %s; line: %s', sequenceId, atomId, line[0] )
                    if errorCount == maxErrorCount+1:
                        ntu.nTerror("And so on")
                    errorCount += 1
                    continue
                #end if
                atmList.append( atm )
            #end for
            if len(atmList) != 2:
                continue
            # Unpack convenience variables.
            atm1 = atmList[0]
            atm2 = atmList[1]
            value = line.float(7)
            if value == 0:
                result().appendPair( (atm1,atm2) )
                continue
            if not value:
                ntu.nTerror("CyanaParser.parseRdcFile: Skipping line without valid value bound on line: [" + line[0]+']')
                continue
            error = line.float(8)
            weight = line.float(9)
            ap = (atm1,atm2)
            r = RDCRestraint( atomPairs= ap, value=value, error=error, weight=weight )
            result.append( r )
            # also store the Candid info if present

        #end for
        if errorCount:
            ntu.nTerror("CyanaParser.parseRdcFile: Found number of errors importing upl file: %s" % errorCount)

        result.rdcFile = rdcFile

        ntu.nTmessage('==> Imported rdc file: %s from %s', result, rdcFile )
        return result
    #end def
    def _fixPdbName(self, name):
        """Fix the digit issue in PDBrecord.name field
        return newName
        """
        newName = moveFirstDigitToEnd(name.strip())
        return newName.strip()
    #end def

    def parsePdbFile(self, pdbFile, silent=False):
        """
        Import PDB records using (modifed) PyMMlib routines.
        Return self or None on success.
        """
        ntu.nTmessage('Parsing pdb file; this may take a while ...')

        self.pdbRecords = PyMMLib.PDBFile(pdbFile)
        if not self.pdbRecords:
            nTerror('CyanaParser.parsePdbFile: parsing PDB-file "%s"', pdbFile)
            return None
        #end if

        # try to match the ATOM and HETATM records, count models
        self.modelCount = 0
        for record in self.pdbRecords:
            recordType = record._name.strip()
            if recordType == 'ATOM' or recordType == 'HETATM':
                chainId    = record.chainID.strip()
                sequenceId = record.resSeq
                residueId  = record.resName.strip()
                atomId     = self._fixPdbName(record.name)

                res = self.findTuple('ResidueRecord', chainId, sequenceId, residueId)
                if res == None and not ( residueId in CYANA_NON_RESIDUES ):        # skip the bloody CYANA non-residue stuff:
                    if not silent:
                        ntu.nTmessage('CyanaParser.parsePdbFile: adding residue (%s, %s, %s)', chainId, sequenceId, residueId)
                    res = ResidueRecord(chainId, sequenceId, residueId)
                    self.appendObject(res)
                #end if

                atm = self.findTuple('AtomRecord', chainId, sequenceId, residueId, atomId)
                if atm == None and not ( residueId in CYANA_NON_RESIDUES ):        # skip the bloody CYANA non-residue stuff:
                    if not silent:
                        ntu.nTmessage('CyanaParser.parsePdbFile: adding atom (%s, %s, %s, %s)', chainId, sequenceId, residueId, atomId)
                    atm = AtomRecord(chainId, sequenceId, residueId, atomId)
                    self.appendObject(atm)
                #end if
                #add a reverse linkage
                record.atom = atm #could be None if pseudo atoms of CYANA_NON_RESIDUES
                record.model = self.modelCount

                if atm != None:
                    atm.pdbRecords.append(record)
                #end if

            elif recordType == 'ENDMDL':
                self.modelCount += 1
            #end if
        #end for

        self.pdbFile = pdbFile
        self.modelCount = max(self.modelCount, 1)
        ntu.nTmessage('==> Parsed %d records (%d models) from %s', len(self.pdbRecords), self.modelCount, pdbFile)
        return self
    #end def

    def _translateTopology(self, residue, topDefList):
        """
        Return translation of topDefList tuples (using CING nomenclature) as a NTlist or None on error
        Adapted from translateTopology in molecule.py to first do a translation to CYANA1/2/3 nomenclature
        """
        if residue == None:
            return None

        result = ntu.NTlist()
        for resdiffIndex,atomName in topDefList:
            # optimized
            res = residue.sibling( resdiffIndex )
            if res == None:
                return None
            if res.typeIdentifier == None:
                return None
            if self.convention == None:
                return None
            # translate the atomName from internal CING nomenclature
            aname = translateAtomName( INTERNAL, res.typeIdentifier.name, atomName, self.convention )
            result.append( ntu.getDeepByKeysOrAttributes( res, aname) )
        #end for
        return result
    #end def

    def _getDihedralAtoms( self, residue, dihedralName ):
        """Return NTlist of AtomRecords for dihdralName of residue or None on error
        """
        # use the getDeepByKeysOrAttributes routine to traverse down the data structure checking for presence at each stage
        topDefList = ntu.getDeepByKeysOrAttributes( residue,'typeIdentifier',dihedralName,'atoms')
        if topDefList == None:
            return None
        return self._translateTopology( residue, topDefList )
    #end def

    def parseAcoFile( self, acoFile ):
        """Read and parse Cyana acoFile
           ( 512 THR   PSI     116.0   148.0)
           convention = CYANA or CYANA2
           return a DihedralRestraintList or None on error
        """
        #maxErrorCount = 50
        #errorCount = 0

        # check for presence of molecule
        if (not self.molecule ):
            ntu.nTerror("parseAcoFile: initialize molecule first")
            return None
        #end if
        ##changed _dir,name,_ext = acoFile.split3() due to error messages stating
        ## str has not attribute split3
        name,ext = acoFile.split('.aco')
        result = DihedralRestraintList( name=name, status='keep')

        for line in AwkLike( acoFile, commentString = '#' , minNF = 5):
          if "type=2" not in line:
          ## filters out CYANA dihedral files that contain multiple dihedrals per residue

            sequenceId   = line.int(1)
            residueId    = line[2].strip()[:3]
            dihedralName = line[3].strip()
            lower        = line.float(4)
            upper        = line.float(5)

            residue = self.findTuple('ResidueRecord', CYANA_CHAIN_ID, sequenceId, residueId)
            if residue == None:
                ntu.nTerror('CyanaParser.parseAcoFile: failed to parse residue %s %s, line %d (%s)',
                             sequenceId, residueId, line.NR, line[0]
                           )
                continue
            #end if

            atoms = self._getDihedralAtoms( residue, dihedralName )
            if atoms == None or (None in atoms):
                ntu.nTerror('CyanaParser.parseAcoFile: failed to decode atoms for dihedral %s, line %d (%s)',
                             dihedralName, line.NR, line[0]
                           )
                continue
            else:
                r = DihedralRestraint( atoms = atoms, lower=lower, upper=upper, angle = dihedralName, residue = residue)
                #print r.format()
                result.append( r )
            #end if
        #end for

        result.fileName = acoFile
        ntu.nTmessage('==> Imported aco file: %s from "%s"', result, acoFile )
        return result
    #end def

    def __str__(self):
        return ntu.sprintf("<CyanaParser: path: %s, residues %d, atoms %d (%s), resonances %d, models %d>",
                           self.rootPath, len(self.sequence), len(self.atoms), self.convention, len(self.resonances),self.modelCount
                           )
    #end def

    def __repr__(self):
        return str(self)
    #end def

#end class
