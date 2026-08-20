"""
Implementation of the CING API's main classes.
Split into 3 for better performance.
"""

from ConfigParser import ConfigParser
from cing import cingPythonCingDir
from cing import cingRoot
from cing import issueListUrl
from cing.Libs.Geometry import violationAngle
from cing.Libs.NTutils import * #@UnusedWildImport
from cing.Libs.disk import copydir
from cing.Libs.disk import remove
from cing.PluginCode.required.reqNih import TALOSPLUS_LIST_STR
from cing.PluginCode.required.reqWhatif import summaryCheckIdList
from cing.STAR.File import File
from cing.core.classes2 import * #@UnusedWildImport
from cing.core.constants import * #@UnusedWildImport
from cing.core.molecule import Atom
from cing.core.molecule import Ensemble
from cing.core.molecule import Molecule
from cing.core.molecule import nTdihedralOpt
from cing.core.molecule import nTdistanceOpt #@UnusedImport
from glob import glob
from glob import glob1
from shutil import rmtree
import tarfile
__version__ = cing.__version__
__date__ = cing.__date__
__author__ = cing.__author__
__copyright__ = cing.__copyright__
__credits__ = cing.__credits__

projects = NTlist()

#: CRV stands for CRiteria Value CRS stands for CRiteria String
CRV_NONE = "-999.9"

#-----------------------------------------------------------------------------
# Cing classes and routines
#-----------------------------------------------------------------------------
# pylint: disable=R0902
#
PeakIndex = 0
class Peak(NTdict, Lister):
    """Peak class:
       Peaks point to resonances
       Resonances point to atoms
       by GV 2007 07 23: added hasHeight, hasVolume attributes to the class
          GV 2007 28 09: Moved from molecule.py to project.py
          GV 19 Jun 08: PeakIndex starts at 0
          GV 19 Jun 08: change height, volume to NTvalue classes
          GV 19 Jun 08: changed hasHeight and hasVolume to methods
    """

    HEIGHT_VOLUME_FORMAT = '%9.2e +- %8.1e'
    HEIGHT_VOLUME_FORMAT2 = '%9.2e'

    def __init__(self,
                  dimension,
                  positions = None,
                  height = NaN, heightError = NaN,
                  volume = NaN, volumeError = NaN,
                  resonances = None,
                  **kwds
                ):

        NTdict.__init__(self, __CLASS__ = 'Peak', **kwds)
        Lister.__init__(self)
#       several external programs need an index
        # Using the global statement pylint: disable=W0603
        global PeakIndex
        self.peakIndex = PeakIndex
        PeakIndex += 1

        self.__FORMAT__ = self.header() + '\n' + \
                           'dimension:  %(dimension)dD\n' + \
                           'positions:  %(positions)s\n' + \
                           'height:     %(height)s\n' + \
                           'volume:     %(volume)s\n' + \
                           'resonances: %(resonances)s\n' + \
                           'rogScore:   %(rogScore)s\n' + \
                           self.footer()

        self.dimension = dimension

        # Copy the positions and resonances argument to assure they become
        # NTlist objects
        if resonances:
            self.resonances = NTlist(*resonances)
        else:
            self.resonances = nTfill(None, dimension)
        #end if

        if positions:
            self.positions = NTlist(*positions)
        else:
            self.positions = nTfill(NaN, dimension)
        #end if

        self.height = NTvalue(height, heightError, Peak.HEIGHT_VOLUME_FORMAT, Peak.HEIGHT_VOLUME_FORMAT2)
        self.volume = NTvalue(volume, volumeError, Peak.HEIGHT_VOLUME_FORMAT, Peak.HEIGHT_VOLUME_FORMAT2)

        self.rogScore = ROGscore()
    #end def

    def decriticize(self):
#        nTdebug("Now in Peak#decriticize")
        self.rogScore.reset()
    #end def

    def isAssigned(self, axis):
        if axis >= self.dimension:
            return False
        if axis >= len(self.resonances):
            return False
        if self.resonances[axis] == None:
            return False
        if self.resonances[axis].atom == None:
            return False
        return True
    #end def

    def getAssignment(self, axis):
        """Return atom instances in case of an assignment or None
        """
        if self.isAssigned(axis):
            return self.resonances[axis].atom
        #end if
        return None
    #end def

    def hasHeight(self):
        return not isNaN(self.height.value)

    def hasVolume(self):
        return not isNaN(self.volume.value)

    def __str__(self):
        #print '>>', self.resonances.zap('atom')
        return sprintf('Peak %4d (%dD)  [%s]   height: %s   volume: %s    Assiged to: %s',
                         self.peakIndex, self.dimension,
                         self.positions.format('%8.2f'),
                         self.height, self.volume,
                         self.resonances.zap('atom').format('%-20s')
                       )
    #end def

    def header(self, mdots = dots):
        """Subclass header to generate using __CLASS__, peakIndex and dots.
        """
        return sprintf('%s %s: %d %s', mdots, self.__CLASS__, self.peakIndex, mdots)
    #end def
#end class


class PeakList(NTlist, ProjectListMember):

    def __init__(self, name, status = 'keep'):
        NTlist.__init__(self)
        ProjectListMember.__init__(self)
        self.name = name
        self.status = status
        self.listIndex = -1 # list is not appended to anything yet
    #end def

    def minMaxDimension(self):
        """Return a tuple of the min and max of the spectral dimensions from all peaks
        usually min == max but not guaranteed.
        Will return (None,None) for empty lists
        """
        minD = None
        maxD = None
        for peak in self:
            minD = min(peak.dimension, minD)
            maxD = max(peak.dimension, maxD)
        return (minD, maxD)

    def peakFromAtoms(self, atoms, onlyAssigned = True):
        """Append a new Peak based on atoms list
           Return Peak instance, or None
        """
        if (None not in atoms):     # None value in atoms indicates atom not present
            if (onlyAssigned and (False in map(Atom.isAssigned, atoms))):
                pass                # Check atom assignments, if only assigned and not all assignments
                                    # present we skip it
            else:                   # all other cases we add a peak
                s = []
                r = []
                for a in atoms:
                    s.append(a.shift())
                    r.append(a.resonances())
                #end for
                peak = Peak(dimension = len(atoms),
                             positions = s,
                             resonances = r
                           )
                self.append(peak)
                return peak
            #end if
        #end if
        return None
    #end def

    def __str__(self):
        return sprintf('<PeakList "%s" (%s,%d)>', self.name, self.status, len(self))
    #end def

    def __repr__(self):
        return str(self)
    #end def

    def format(self):  # pylint: disable=W0221
        s = sprintf('%s PeakList "%s" (%s,%d,%s) %s\n', dots, self.name, self.status, len(self), self.rogScore, dots)
        for peak in self:
            s = s + str(peak) + '\n'
        #end for
        return s
    #end def

    def save(self, path = None):
        """
        Create a SML file
        Return self or None on error
        """
        if not path:
            path = self.objectPath
        if self.SMLhandler.toFile(self, path) != self:
            nTerror('PeakList.save: failed creating "%s"', path)
            return None
        #end if

        nTdetail('==> Saved %s to "%s"', self, path)
        return self
    #end def

    def rename(self, newName):
        return self.projectList.rename(self.name, newName)
    #end def
#end class


def getAtomsFromAtomPairs(atomPairs):
    result = []
    for atomPair in atomPairs:
        for atom in atomPair:
            for real_atm in atom.realAtoms():
                result.append(real_atm)
    return result

class Restraint(NTdict):
    """
    Super class for DistanceRestraint etc.
    On initialization the atom pairs will be tested for validity and reported in self.isValid
    """
    def __init__(self, lower, upper, **kwds):
        NTdict.__init__(self, lower = lower,
                              upper = upper,
                              **kwds
                       )
        self.__CLASS__ = None     # Set by sub class.
        self.id = -1       # Undefined index number
        self.violations = None     # list with violations for each model, None indicates no analysis done
        self.violCount1 = 0        # Number of violations over 1 degree (0.1A)
        self.violCount3 = 0        # Number of violations over 3 degrees (0.3A)
        self.violCount5 = 0        # Number of violations over 5 degrees (0.5A)
        self.violMax = 0.0      # Maximum violation
        self.violAv = 0.0      # Average violation
        self.violSd = 0.0      # Sd of violations
        self.isValid = True
        self.rogScore = ROGscore()
    #end def
    def __str__(self):
        return '<%s %d>' % (self.__CLASS__, self.id)
    #end def

    def decriticize(self):
#        nTdebug("Now in Restraint#%s" % getCallerName())
        self.rogScore.reset()
    #end def

    def getModelCount(self):
        """Iterate over the atoms until an atom is found that returns not a None for getModelCount.
        Return 0 if it doesn't have any models or None on error.
        """
        modelCount = None

        if self.__CLASS__ == DR_LEVEL or self.__CLASS__ == RDC_LEVEL:
#            lAtom  = len(self.atomPairs)
            for atompair in self.atomPairs:
                for atom in atompair:
                    if not isinstance(atom, Atom):
                        nTerror("Failed to get atom in atom pair in %s for %s" % (getCallerName(), self))
                        return None
                    modelCount = atom.getModelCount()
                    if modelCount != None:
                        return modelCount
        else:
#            lAtom  = len(self.atoms)
            for atom in self.atoms:
                modelCount = atom.getModelCount()
                if modelCount != None:
                    return modelCount
#        nTwarning("%s.getModelCount returned None for all %d atom(pair)s; giving up." % (self.__CLASS__, lAtom))
        return None
    #end def

    def isValidForAquaExport(self):
        """Determine if the restraint can be exported to Aqua."""
        nTerror("Restraint.isValidForAquaExport needs to be overriden.")
    #end def

    def listViolatingModels(self, cutoff = 0.3):
        """
        Examine for violations larger then cutoff, return list of violating models or None on error
        Requires violations attribute (obtained with calculateAverage method.
        """
        if not self.has_key('violations'):
            return None

        violatingModels = NTlist()
        for i in range(len(self.violations)):
            if (math.fabs(self.violations[i]) > cutoff):
                violatingModels.append(i)
            #end if
        #end for

        return violatingModels
    #end def
# end class Restraint


class DistanceRestraint(Restraint):
    """DistanceRestraint class:
       atomPairs: list of (atom_1,atom_2) tuples,
       lower and upper bounds.
    """
    STATUS_SIMPLIFIED = 'simplified'
    STATUS_NOT_SIMPLIFIED = 'not simplified'
    STATUS_DEASSIGNED = 'deassigned'
    STATUS_NOT_DEASSIGNED = 'not deassigned'
    STATUS_REMOVED_DUPLICATE = 'removed duplicate'
    STATUS_NOT_REMOVED_DUPLICATE = 'not removed duplicate'
    # The maximum number of atom pairs expected before it will be treated normally.
    # This is to prevent HADDOCK AIR restraints to slow CING to a crawl as per issue 324.
    MAX_ATOM_PAIRS_EXPECTED = 1000

#    def __init__( self, atomPairs=[], lower=0.0, upper=0.0, **kwds ):
    def __init__(self, atomPairs = NTlist(), lower = None, upper = None, **kwds):
        Restraint.__init__(self, lower = lower, upper = upper, **kwds)
        self.__CLASS__ = DR_LEVEL
        self.atomPairs = NTlist()
        self.distances = None     # list with distances for each model; None: not yet defined
        self.av = None      # Average distance
        self.sd = None      # sd on distance
        self.min = None      # Minimum distance
        self.max = None      # Max distance
        self.violCountLower = 0    # Lower-bound violations
        self.violUpperMax = 0.0    # Max violation over upper bound
        self.violLowerMax = 0.0    # Max violation over lower bound

        self.duplicates = NTlist() # NTlist instance with DistanceRestraint instances considered duplicates; TODO: check this code
        self.error = False    # Indicates if an error was encountered when analyzing restraint

        for pair in atomPairs:
            if self.appendPair(pair):
#                nTdebug('resetting self.isValid')
                self.isValid = False
                return
        #end for
    #end def

    def isValidForAquaExport(self):
        """Determine if the restraint can be exported to Aqua.
        Simplified to checking if all partners have at least one assigned atom.

        E.g. of a valid self.atomPairs
        [ [ HA ] [ HB,HC ],
          [ HA ] [ HD ] ]
        """
        if not self.atomPairs:
#            nTdebug("Failed to find any atom pair in %s" % self)
            return False
        for _i, atomPair in enumerate(self.atomPairs):
            if not atomPair: # eg [ HA ] [ HB,HC ]
#                nTdebug("Failed to find any atomList (should always be 2 present) in atompair %d of:\n%s" % (i,self))
                return False
            for _j, atomList in enumerate(atomPair):
                if not atomList: # eg [ HB,HC ]
#                    nTdebug("Failed to find any atom in atomList (%d,%d) of %s" % (i,j,self))
                    return False
                for _k, atom in enumerate(atomList):
                    if not atom: # eg HB
#                        nTdebug("Failed to find atom in atomList (%d,%d,%d) of %s" % (i,j,k,self))
                        return False
        return True
    #end def


    def criticize(self, project):
        """Only the self violations,violMax and violSd needs to be set before calling this routine"""
        self.rogScore.reset()
#        nTdebug( '%s' % self )
        if (project.valSets.DR_MAXALL_BAD != None) and (self.violMax >= project.valSets.DR_MAXALL_BAD):
            comment = 'violMax: %7.2f' % self.violMax
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_RED, comment)
        elif (project.valSets.DR_MAXALL_POOR != None) and (self.violMax >= project.valSets.DR_MAXALL_POOR):
            comment = 'violMax: %7.2f' % self.violMax
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_ORANGE, comment)
        if (project.valSets.DR_THRESHOLD_OVER_POOR != None) and (project.valSets.DR_THRESHOLD_FRAC_POOR != None):
            fractionAbove = getFractionAbove(self.violations, project.valSets.DR_THRESHOLD_OVER_POOR)
            if fractionAbove >= project.valSets.DR_THRESHOLD_FRAC_POOR:
                comment = 'fractionAbove: %7.2f' % fractionAbove
    #            nTdebug(comment)
                self.rogScore.setMaxColor(COLOR_ORANGE, comment)
        if (project.valSets.DR_THRESHOLD_OVER_BAD != None) and (project.valSets.DR_THRESHOLD_FRAC_BAD != None):
            fractionAbove = getFractionAbove(self.violations, project.valSets.DR_THRESHOLD_OVER_BAD)
            if fractionAbove >= project.valSets.DR_THRESHOLD_FRAC_BAD:
                comment = 'fractionAbove: %7.2f' % fractionAbove
    #            nTdebug(comment)
                self.rogScore.setMaxColor(COLOR_RED, comment)
        if (project.valSets.DR_RMSALL_BAD != None) and (self.violSd >= project.valSets.DR_RMSALL_BAD):
            comment = 'violSd: %7.2f' % self.violSd
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_RED, comment)
        elif (project.valSets.DR_RMSALL_POOR != None) and (self.violSd >= project.valSets.DR_RMSALL_POOR):
            comment = 'violSd: %7.2f' % self.violSd
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_ORANGE, comment)



        if project.valSets.FLAG_MISSING_COOR:
            #modelCount = self.getModelCount()
            for atm1, atm2 in self.atomPairs:
                atms1 = atm1.realAtoms()
                atms2 = atm2.realAtoms()
                for a in atms1 + atms2:
                    if a and a.hasMissingCoordinates(): # gv has mase this into a method because the getModelCount()
                    # can crash when reading back NRG dataset because of their
                    # incompleteness
                    #if len(a.coordinates) < modelCount:
                        msg = "Missing coordinates (%s)" % a.toString()
#                        nTdebug(msg)
                        self.rogScore.setMaxColor(COLOR_RED, msg)
                    #end if
                #end for
            #end for
        # end if
    #end def

    def simplify(self):
        """
        Return True on error.

        Routine is iterative itself because both sides may contain ambis to collapse and remove.
        """
        atomPairCount = len(self.atomPairs)
        if atomPairCount > DistanceRestraint.MAX_ATOM_PAIRS_EXPECTED: # Happens for entry 2bgf as per issue 324.
#            nTdebug('In %s; skipping restraint %s with %s atom pairs which is more than the maximum expected: %s' % (
#                getCallerName(), self, atomPairCount, DistanceRestraint.MAX_ATOM_PAIRS_EXPECTED))
            return self.STATUS_NOT_SIMPLIFIED
        # end if

        statusOverall = self.STATUS_NOT_SIMPLIFIED
        status = self.STATUS_SIMPLIFIED
        while status == self.STATUS_SIMPLIFIED:
            status = self.simplifyForFc()
            if status == self.STATUS_SIMPLIFIED:
                statusOverall = status
#                nTdebug("simplified restraint %s" % self)
            elif status == self.STATUS_NOT_SIMPLIFIED:
                pass
#                nTdebug("not simplified restraint %s" % self)
            else:
                nTerror("Encountered an error simplifying restraint %s" % self)
                return True
            # end if
        # end while

        if self.removeDuplicateAtomPairs():
            nTerror("Encountered an error in removeDuplicateAtomPairs restraint %s" % self)
            return True
        # end if
        while self.removeDuplicateAtomPairs2() == self.STATUS_REMOVED_DUPLICATE:
            pass
#            nTdebug("Removed duplicate")
        # end if
        return statusOverall
    #end def


    def deassignStereospecificity(self):
        """If the restraint involves a stereo specifically assignable atom then expand the list to include all
        getStereoPartner's. Of course if the restraint is between partners then the restraint
        becomes useless but will be generated. E.g. Gly HA2 to HA3 will become Gly QA to QA.

        LEU MD1 -> QD
        PHE QD  -> QR (which is a problem for the xplor conversion as of yet.)
        Return None on error.
        STATUS_DEASSIGNED = 'deassigned'
        STATUS_NOT_DEASSIGNED = 'not deassigned'
        """

        nTdebug('Starting 123 deassignStereospecificity for %s' % ( self ) )
        isDeassigned = False

        for atomPairIdx in range(len(self.atomPairs)):
            atomPair = self.atomPairs[atomPairIdx]
            for atomIdx in range(2):
                atomOrg = atomPair[atomIdx]
                pseudoAtom = atomOrg.pseudoAtom()
                if not pseudoAtom:
                    continue
                # Deal with immutable tuples.
                if atomIdx == 0:
                    newTuple = ( pseudoAtom, atomPair[1])
                else:
                    newTuple = ( atomPair[0], pseudoAtom )
                atomPair = newTuple
                self.atomPairs[atomPairIdx] = atomPair
                isDeassigned = True
                nTdebug('Replaced %s by %s' % ( atomOrg, atomPair[atomIdx] ) )
        if isDeassigned:
            return self.STATUS_DEASSIGNED
        return self.STATUS_NOT_DEASSIGNED

    def simplifyForFc(self):
        """FC likes to split Val QQG in QG1 and 2 making it appear to be an ambiguous OR typed XPLOR restraint
        were it is not really one. Undone here.

        Returns:
        None                     error.
        STATUS_NOT_SIMPLIFIED    no simplifications done
        STATUS_SIMPLIFIED        simplifications done

        In the code:

        j stands for the index of the atomPair of the outer loop that might be removed upon simplification.
        i stands for the index of the atomPair of the inner loop that is compared to and that might be modified to include atoms from atomPair j.
        """
#        nTdebug('Starting simplifyForFc for\n:%r' % ( self ) )
        atomPairIdxJ = len(self.atomPairs) # starting from the end.
        while atomPairIdxJ > 1:
            atomPairIdxJ -= 1
            atomPairJ = self.atomPairs[atomPairIdxJ]
            atomPairJset = set(atomPairJ) # Important to use api of unsorted atoms in pair (left right will not matter)
            atom0J = atomPairJ[0]
            atom1J = atomPairJ[1]

#            nTdebug('For atomPairIdxJ %d using atoms J %s and %s' % ( atomPairIdxJ, atom0J, atom1J) )
            # speed up check on J as an early abort clause.
            if not (atom0J.hasPseudoAtom() or atom1J.hasPseudoAtom()):
                if not (atom0J.getPseudoOfPseudoAtom() or atom1J.getPseudoOfPseudoAtom()):
#                    nTdebug('Skipping restraint without pseudo representing J atoms')
                    continue

            for atomPairIdxI in range(atomPairIdxJ): # Compare only with the previous atom pairs
                atomPairI = self.atomPairs[atomPairIdxI]
                _atom0I = atomPairI[0]
                _atom1I = atomPairI[1]
#                nTdebug('    Using atoms I %s and %s' % ( atom0I, atom1I) )
                atomPairIset = set(atomPairI)
                atomPairIntersection = atomPairIset.intersection(atomPairJset)
                if not atomPairIntersection:
#                    nTdebug('    No intersection')
                    continue

#                 At this point it is certain that there is an intersection of at least one atom between the two pairs.
                if len(atomPairIntersection) != 1:
#                    nTdebug('More than one atom in atom set intersection: %s' % atomPairIntersection)
                    continue

                atomInCommon = atomPairIntersection.pop() # get arbitrary element of set.
                atomIinCommonIdx = 0
                atomJinCommonIdx = 0
                atomItoMergeIdx = 1
                atomJtoMergeIdx = 1
                if atomPairI[atomIinCommonIdx] != atomInCommon:
                    atomIinCommonIdx = 1
                    atomItoMergeIdx = 0
                if atomPairJ[atomJinCommonIdx] != atomInCommon:
                    atomJinCommonIdx = 1
                    atomJtoMergeIdx = 0

                # Now we know which atoms are in common and consequently the others should be tried to merge.
#                nTdebug('    atominCommonIdx I %d and J %d for %s' % ( atomIinCommonIdx, atomJinCommonIdx, atomInCommon) )

                atomItoMerge = atomPairI[atomItoMergeIdx]
                atomJtoMerge = atomPairJ[atomJtoMergeIdx]

                atomIinCommon = atomPairI[atomIinCommonIdx]
                atomJinCommon = atomPairJ[atomJinCommonIdx]

#                nTdebug('    atomIinCommon %s == atomJinCommon %s' % ( atomIinCommon, atomJinCommon ))
                if atomIinCommon != atomJinCommon:
                    nTcodeerror('    atoms toMerge I %s and J %s differ.' % ( atomItoMerge, atomJtoMerge) )
                    continue
                # end if

                if atomItoMerge.getStereoPartner() != atomJtoMerge:
#                    nTdebug('    atoms toMerge I %s and J %s have different parent if at all related.' % ( atomItoMerge, atomJtoMerge) )
                    continue
                # end if

                pseudoOfAtom = atomItoMerge.pseudoAtom()
                if not pseudoOfAtom:
#                    nTdebug('    no pseudo for this atom %s' % atomItoMerge)
                    pseudoOfAtom = atomItoMerge.getPseudoOfPseudoAtom()
                    if not pseudoOfAtom:
                        nTwarning('    no pseudo of pseudoatom %s' % atomItoMerge) # happens in 1y0j for <Atom A.VAL205.CG1>
                        continue
                    # end if
                # end if

#                nTdebug( "    New pop atom: %s" % pseudoOfAtom)
                # Change I maintaining order
                atomPairINewList = list(atomPairI)
                atomPairINewList[atomItoMergeIdx] = pseudoOfAtom
                self.atomPairs[atomPairIdxI] = tuple(atomPairINewList)
#                nTdebug("Now self.atomPairs[atomPairIdxI]: %s" % str(self.atomPairs[atomPairIdxI]))
                # Remove J
#                nTdebug("Removing self.atomPairs[atomPairIdxJ]: %s" % str(self.atomPairs[atomPairIdxJ]))
                del self.atomPairs[atomPairIdxJ]
                # Return quickly to keep code to the left (keep it simple).
#                nTdebug('Simplified.')
                return self.STATUS_SIMPLIFIED
            # end for
        # end while
#        nTdebug('Not simplified.')
        return self.STATUS_NOT_SIMPLIFIED
    # end def

    def removeDuplicateAtomPairs(self):
        """
        Used in simplify.

        Returns:
        True                     error.

        In the code:

        j stands for the index of the atomPair of the outer loop that might be removed upon removal.
        i stands for the index of the atomPair of the inner loop that is compared to.
        """

#        nTdebug('Starting %s for %s' % ( getCallerName(), self ) )
        atomPairIdxJ = len(self.atomPairs) # starting from the end.
        while atomPairIdxJ > 1:
            atomPairIdxJ -= 1
            atomPairJ = self.atomPairs[atomPairIdxJ]
            atomPairJset = set(atomPairJ) # Important to use api of unsorted atoms in pair (left right will not matter)
            _atom0J = atomPairJ[0]
            _atom1J = atomPairJ[1]

#            nTdebug('For atomPairIdxJ %d using atoms J %s and %s' % ( atomPairIdxJ, atom0J, atom1J) )

            for atomPairIdxI in range(atomPairIdxJ): # Compare only with the previous atom pairs
                atomPairI = self.atomPairs[atomPairIdxI]
                _atom0I = atomPairI[0]
                _atom1I = atomPairI[1]
#                nTdebug('    Using atoms I %s and %s' % ( atom0I, atom1I) )
                atomPairIset = set(atomPairI)
                atomPairIntersection = atomPairIset.intersection(atomPairJset)
                if not atomPairIntersection:
#                    nTdebug('    No intersection')
                    continue
                if len(atomPairIntersection) != 2:
#                    nTdebug('Only one atom in atom set intersection: %s' % atomPairIntersection)
                    continue
#                nTdebug("Removing self.atomPairs[atomPairIdxJ]: %s" % str(self.atomPairs[atomPairIdxJ]))
                del self.atomPairs[atomPairIdxJ]
            # end for
        # end while
        return
    # end def

    def removeDuplicateAtomPairs2(self):
        """
        Used in simplify.

        This code is more advanced than the above removeDuplicateAtomPairs2 in that it will also
        check when pseudos are contained in other pseudos. The widest will be retained.
        E.g.
        For 1a24
        783.00    A    20    PRO    QB    A    23    LEU    MD1   3.20    7.90    2.96    0.56    2.56    3.35    0.32    0.45    0.64    0    0    0
        783.01    A    20    PRO    QB    A    23    LEU    QD    3.20    7.90    2.96    0.56    2.56    3.35    0.32    0.45    0.64    0    0    0
        will be truncated to:
        For 1a24
        783       A    20    PRO    QB    A    23    LEU    QD    3.20    7.90    2.96    0.56    2.56    3.35    0.32    0.45    0.64    0    0    0

        Watch for e.g. intraresidual LEU with multiple atompairs:
        HB2 QD
        QB  MD1
        This can not be truncated whereas:
        QB  MD1
        HB2 QD
        HB3 QD
        can be to simply QB QD.
        Therefore previous routines should already have cleaned up to:
        QB  MD1
        QB  QD
        so that this routine will do the final collapse to QB QD.

        The ordering is irrelevant but must always be maintained.


        pseudo code:
        loop over atompairs i,j
                atomset0i
                if(
                   ( atomset0i.issuperset(atomset0j) and atomset1i.issuperset(atomset1j)) or
                   ( atomset0i.issuperset(atomset1j) and atomset1i.issuperset(atomset0j))    )
                    remove j and return
                elif(
                   ( atomset0j.issuperset(atomset0i) and atomset1j.issuperset(atomset1i)) or
                   ( atomset0j.issuperset(atomset1i) and atomset1j.issuperset(atomset0i))    )
                    remove i and return

        Returns:
        True                     error.
        STATUS_REMOVED_DUPLICATE = 'removed duplicate'
        STATUS_NOT_REMOVED_DUPLICATE = 'not removed duplicate'

        In the code:

        j stands for the index of the atomPair of the outer loop that might be removed upon removal.
        i stands for the index of the atomPair of the inner loop that is compared to.
        """

#        nTdebug('Starting %s for %s' % ( getCallerName(), self ) )

        n = len(self.atomPairs)
        for atomPairIdxJ in range(n-1):
            atomPairJ = self.atomPairs[atomPairIdxJ]
#            atomPairJset = set(atomPairJ) # Important to use api of unsorted atoms in pair (left right will not matter)
            atom0J = atomPairJ[0]
            atom1J = atomPairJ[1]
            atomset0J = set( atom0J.realAtoms() )
            atomset1J = set( atom1J.realAtoms() )

#            nTdebug('For atomPairIdxJ %d using atoms J %s and %s' % ( atomPairIdxJ, atom0J, atom1J) )

            for atomPairIdxI in range(atomPairIdxJ+1,n): # Compare only with the next atom pairs
                atomPairI = self.atomPairs[atomPairIdxI]
                atom0I = atomPairI[0] #@UnusedVariable
                atom1I = atomPairI[1] #@UnusedVariable
#                nTdebug('    Using atoms I %s and %s' % ( atom0I, atom1I) )

                atomset0I = set( atom0I.realAtoms() )
                atomset1I = set( atom1I.realAtoms() )
                if(
                   ( atomset0I.issuperset(atomset0J) and atomset1I.issuperset(atomset1J)) or
                   ( atomset0I.issuperset(atomset1J) and atomset1I.issuperset(atomset0J))    ):
#                    nTdebug("Removing self.atomPairs[atomPairIdxJ]: %s" % str(self.atomPairs[atomPairIdxJ]))
                    del self.atomPairs[ atomPairIdxJ ]
                    return self.STATUS_REMOVED_DUPLICATE
                elif(
                   ( atomset0J.issuperset(atomset0I) and atomset1J.issuperset(atomset1I)) or
                   ( atomset0J.issuperset(atomset1I) and atomset1J.issuperset(atomset0I))    ):
#                    nTdebug("Removing self.atomPairs[atomPairIdxI]: %s" % str(self.atomPairs[atomPairIdxI]))
                    del self.atomPairs[ atomPairIdxI ]
                    return self.STATUS_REMOVED_DUPLICATE
                # end if
            # end for
        # end while
        return self.STATUS_NOT_REMOVED_DUPLICATE
    # end def

    def appendPair(self, pair):
        """ pair is a (atom1,atom2) tuple

        check if atom1 already present, keep order
        otherwise: keep atom with lower residue index first
        Return True on error.
        """
        # GV says; order needs to stay: is being used for easier
        # (manual) analysis.


        if pair[0] == None or pair[1] == None:
            nTerror('DistanceRestraint.appendPair: invalid pair %s', str(pair))
            return True
        #end if

#        missesId = False
        for atom in pair:
            if not hasattr(atom, 'id'): # happens for 1f8h and LdCof (imported from CYANA data).
#                nTwarning('DistanceRestraint.appendPair: invalid pair %s for atom missing id: %s' % (str(pair), str(atom)))
#                missesId = True
                return True
        #end if

#        if missesId:
#            self.atomPairs.append((pair[0], pair[1]))
#        else:
        # gv 24 Jul: just use atoms id, they are unique and ordered
        if pair[0].id < pair[1].id:
            self.atomPairs.append((pair[0], pair[1]))
        else:
            self.atomPairs.append((pair[1], pair[0]))
    #end def

    def classify(self):
        """
        Return 0,1,2,3 depending on sequential, intra-residual, medium-range or long-range
        Simply ignore ambiguous assigned NOEs for now and take it as the first atom pair

        return -1 on error
        """
        if not self.atomPairs or len(self.atomPairs) < 1:
            return - 1

        atm1, atm2 = self.atomPairs[0]
        if atm1 == None or atm2 == None:
            return - 1

        if atm1.residue.chain != atm2.residue.chain:
            return 3
        elif atm1.residue == atm2.residue:
            return 0
        else:
            r1 = atm1.residue
            r2 = atm2.residue
            if not (r1.chain and r2.chain):
                # Deals with removed residues which don't have a parent anymore.
                return -1
            delta = int(math.fabs(r1.chain._children.index(r1) - r2.chain._children.index(r2)))
            if delta == 1:
                return 1
            elif delta > 4:
                return 3
            else:
                return 2
            #end if
        #end if
    #end def

    def isAmbiguous(self):
        return len(self.atomPairs) > 1
    #end def

    def calculateAverage(self):
        """
        Calculate R-6 average distance and violation
        for each model.
        return (av,sd,min,max) tuple, or (None, None, None, None) on error
        Important to set the violations to 0.0 if no violations were found.
        In that case the s.d. may remain None to indicate undefined.
        """

#        nTdebug('calculateAverage: %s' % self)
        self.error = False    # Indicates if an error was encountered when analyzing restraint

        modelCount = self.getModelCount()
        if not modelCount:
#            nTdebug('DistanceRestraint.calculateAverage: No structure models (%s)', self)
            self.error = True
            return (None, None, None, None)
        #end if

        self.distances = NTlist() # list with distances for each model
        self.av = None      # Average distance
        self.sd = None      # sd on distance
        self.min = None      # Minimum distance
        self.max = None      # Max distance
        self.violations = NTlist() # list with violations for each model INCLUDING non violating models!
        self.violCount1 = 0        # Number of violations > 0.1A
        self.violCount3 = 0        # Number of violations > 0.3A
        self.violCount5 = 0        # Number of violations > 0.5A
        self.violCountLower = 0    # Number of lower-bound violations over 0.1A
        self.violMax = 0.0         # Maximum violation
        self.violUpperMax = 0.0    # Max violation over upper bound
        self.violLowerMax = 0.0    # Max violation over lower bound
        self.violAv = 0.0          # Average violation
        self.violSd = None         # Sd of violations
        self.violSum = 0.0         # Sum of violations
        self.distances = nTfill(0.0, modelCount) #list with actual effective distances

        models = range(modelCount)
        i = 0
        atm1, atm2 = None, None # helping pylint.
        for atm1, atm2 in self.atomPairs:

            # GV says: Check are done to prevent crashes upon rereading
            # datasets with floating/adhoc residues/atoms

            # skip trivial cases
            if atm1 == atm2:
                continue
            if atm1 == None or atm2 == None:
                continue
            #expand pseudoatoms
            atms1 = atm1.realAtoms()
            if atms1 == None:
                #nTdebug('DistanceRestraint.calculateAverage: %s.realAtoms() None (%s)', atm1, self)
                continue
            atms2 = atm2.realAtoms()
            if atms2 == None:
                #nTdebug('DistanceRestraint.calculateAverage: %s.realAtoms() None (%s)', atm2, self)
                continue
            for a1 in atms1:
                #print '>>>', a1.format()
                if len(a1.coordinates) == modelCount:
                    for a2 in atms2:
                        #print '>>', atm1, a1, atm2, a2
                        i = 0
                        if len(a2.coordinates) == modelCount:
                            for i in models:
                                self.distances[i] += Rm6dist(a1.coordinates[i].e, a2.coordinates[i].e)
                            #end for
                        else:
#                            self.distances[0] = 0.0
                            i = 0
                            self.error = True
                            break
                        #end if
                    #end for
                else:
#                    self.distances[0] = 0.0
                    i = 0
                    self.error = True
                    break
                #end if
            #end for
        #end for
        if self.error:
            msg = "AtomPair (%s,%s) model %d without coordinates" % (atm1.toString(), atm2.toString(), i)
#            nTdebug(msg)
            self.rogScore.setMaxColor(COLOR_RED, msg)
            return (None, None, None, None)
        #end if

        # Calculate R-6 distances
        for i in models:
            if self.distances[i] > 0.0:
                self.distances[i] = math.pow(self.distances[i], -0.166666666666666667)
            #end if
        #end for

        self.av, self.sd, self.n = nTaverage(self.distances)
        self.min = min(self.distances)
        self.max = max(self.distances)

        # calculate violations
        for d in self.distances:
            if d < self.lower:
                self.violations.append(d - self.lower)
            elif (d > self.upper):
                if self.upper == None: # Happens for entry 1but
                    self.violations.append(0.0)
                else:
                    self.violations.append(d - self.upper)
            else:
                self.violations.append(0.0)
            #end if
        #end for

        # analyze violations
        for v in self.violations:
            if (v > 0.5):
                self.violCount5 += 1
            if (v > 0.3):
                self.violCount3 += 1
            if (v > 0.1):
                self.violCount1 += 1
            if (v < -0.1):
                self.violCountLower += 1
        #end for
        if self.violations:
            vAbs = map(math.fabs, self.violations)
            self.violAv, self.violSd, _n = nTaverage(vAbs)
            self.violMax = max(vAbs)
            self.violSum = sum(vAbs)
            self.violUpperMax = max(self.violations)
            self.violLowerMax = min(self.violations)
            if self.violLowerMax < 0.0:
                self.violLowerMax = math.fabs(self.violLowerMax)
            else:
                self.violLowerMax = 0.0
        #end if

        return (self.av, self.sd, self.min, self.max)
    #end def

    def _names(self):
        """
        Internal routine: generate string from atomPairs
        """
        s = ''
        for p in self.atomPairs:
#            s = s + sprintf('(%-11s - %11s)   ', p[0].cName(1), p[1].cName(1))
            s = s + sprintf('(%-13s - %13s)   ', p[0].cName(2), p[1].cName(2)) # include chain id.
        #end for
        return s.strip()
    #end def

    def format(self, allowHtml = False):
        msg = sprintf('%-25s %-6s (Target: %s %s)  (Models: min %s  av %s+-%s  max %s) ' + \
                    '(Violations: av %4.2f max %4.2f counts l,0.1,0.3,0.5:%2d,%2d,%2d,%2d) %s',
                     str(self), self.rogScore,
                     val2Str(self.lower, "%4.1f", 4),
                     val2Str(self.upper, "%4.1f", 4),
                     val2Str(self.min, "%4.1f", 4),
                     val2Str(self.av, "%4.2f", 4),
                     val2Str(self.sd, "%4.1f", 4),
                     val2Str(self.max, "%4.1f", 4),
                     self.violAv, self.violMax,
                     self.violCountLower, self.violCount1, self.violCount3, self.violCount5,
                     self._names()
                    )
        if allowHtml:
            msg = addPreTagLines(msg)
        return msg
    #end def
#end class

# Too many ancestors (8/7) pylint: disable=R0901
class DistanceRestraintList(RestraintList):
    """
    Class based on NTlist that holds distanceRestraints.
    Also manages the "id's".
    Sort  by item of DistanceRestraint Class
    """
    # use the same spelling through out.
    def __init__(self, name, status = 'keep'):
        RestraintList.__init__(self, name, status = status)
        self.__CLASS__ = DRL_LEVEL
        self.hBond = False       # hBond: fix to keep information about hBond restraints from CCPN
        self.violCountLower = 0   # Total lower-bound violations over 0.1 A
        self.violCount1 = 0
        self.violCount3 = 0
        self.violCount5 = 0
        self.violMax = 0.0   # Maximum violation
        self.violUpperMax = 0.0   # Max violation over upper bound
        self.violLowerMax = 0.0   # Max violation over lower bound

        # partitioning in intra, sequential, medium-range and long-range, ambigous
        self.intraResidual = NTlist()
        self.sequential = NTlist()
        self.mediumRange = NTlist()
        self.longRange = NTlist()
        self.ambiguous = NTlist()

        # Duplicate analysis
        self.uniqueDistancesCount = 0       # count of all defined distance restraints
        self.withoutDuplicates = NTlist()   # list of all restraints without duplicates
        self.withDuplicates = NTlist()      # list of all restraints with duplicates
    #end def

    def criticize(self, project, toFile = True):
        """
        Criticize restraints of this list; infer own ROG score from individual restraints.
        """
#        nTdebug('DistanceRestraintList.criticize %s', self)
        self.rogScore.reset()

        for dr in self:
            dr.criticize(project)
            self.rogScore.setMaxColor(dr.rogScore.colorLabel, comment = 'Cascaded from: %s' % dr)
        if toFile:
            path = project.moleculePath('analysis', self.name + '.txt')
            f = file(path, 'w')
            fprintf(f, '%s\n\n', self.format())
            for dr in self:
                fprintf(f, '%s\n', dr.format())
            #end for
            f.close()
            nTdetail('==> Analyzing %s, output to %s', self, path)
        #end if
    #end def

    def simplify(self):
        """Look at Wattos code for a full set of code that does any simplification"""
        for dr in self:
            dr.simplify()

    def deassignStereospecificity(self, max_messages=100):
        """Deassign all"""
        msgHol = MsgHoL()
        for dr in self:
            strDrAtomPairs = str(dr.atomPairs)
            status = dr.deassignStereospecificity()
            if status == DistanceRestraint.STATUS_DEASSIGNED:
                msgHol.appendMessage("%s deassigned to:%s" % (strDrAtomPairs,str(dr.atomPairs)))
        msgHol.showMessage(max_messages=max_messages)

    def analyze(self):
        """
        Calculate averages for every restraint.
        Partition restraints into classes.
        Analyze for duplicate restraints.

        Return <rmsd>, sd and total violations over 0.1, 0.3, 0.5 A as tuple
        or (None, None, None, None, None) on error
        """

#        nTdebug('DistanceRestraintList.analyze: %s', self)

        if (len(self) == 0):
            # happens for entry 2k0e imported from CCPN. Has unlinked restraints.
            nTdebug('DistanceRestraintList.analyze: "%s" empty list'% self.name )
            return (None, None, None, None, None)
        #end if

        modelCount = self.getModelCount()
        if not modelCount:
            nTerror('DistanceRestraintList.analyze: "%s" modelCount %s' % (self.name, modelCount))
            return (None, None, None, None, None)
        #end if

        # check for duplicate
        self.findDuplicates()

        self.rmsd = nTfill(0.0, modelCount)
        self.violCount1 = 0
        self.violCount3 = 0
        self.violCount5 = 0
        count = 0
        self.errors = NTlist() # Store reference to restraints with calc problems

        # partitioning in intra, sequential, medium-range and long-range, ambigous
        self.intraResidual = NTlist()
        self.sequential = NTlist()
        self.mediumRange = NTlist()
        self.longRange = NTlist()
        self.ambiguous = NTlist()

        for dr in self:
            dr.calculateAverage()
            if dr.error:
                self.errors.append(dr)
            else:
                self.violCountLower += dr.violCountLower
                self.violCount1 += dr.violCount1
                self.violCount3 += dr.violCount3
                self.violCount5 += dr.violCount5
                for i in range(0, modelCount):
                    if dr.violations[i]:
                        self.rmsd[i] += dr.violations[i] * dr.violations[i]
                #end for
                count += 1
            #end if

            c = dr.classify()
            if c == 0:
                self.intraResidual.append(dr)
            elif c == 1:
                self.sequential.append(dr)
            elif c == 2:
                self.mediumRange.append(dr)
            elif c == 3:
                self.longRange.append(dr)
            #end if
            if dr.isAmbiguous():
                self.ambiguous.append(dr)
        #end for

        #Set max violations
        for p in ['violMax', 'violUpperMax', 'violLowerMax']:
            myList = self.zap(p)
            myList.sort()
            setattr(self, p, myList[-1])
        #end for

        for i in range(0, modelCount):
            if count:
                if self.rmsd[i]:
                    self.rmsd[i] = math.sqrt(self.rmsd[i] / count)
                else:
                    self.rmsd[i] = None
                #end if
            #end if
        #end for
        self.rmsdAv, self.rmsdSd, _n = nTaverage(self.rmsd)
        return (self.rmsdAv, self.rmsdSd, self.violCount1, self.violCount3, self.violCount5)
    #end def

    def findDuplicates(self):
        """
        Find the duplicate entries in the list.
        TODO: check code with unit check.
        Checked by hand for entry 1cjg on 2012-04-23
        """
        pairs = {} # Hashed by a sorted tuple of atom pairs. The value is a list of drs.
        for dr in self:
            dr.atomPairs.sort() # improves matching for ambiguous restraints
            t = tuple(dr.atomPairs)
            # Reminder: setdefault only creates the value [] if t doesn't exist yet.
            pairs.setdefault(t, [])
            pairs[t].append(dr)
        #end for
        self.uniqueDistancesCount = len(pairs) # Number of different atom pairs.
#        nTdebug("Found %d unique tuples of atom pairs (uniqueDistancesCount)" % self.uniqueDistancesCount)
        self.withoutDuplicates = NTlist()
        self.withDuplicates = NTlist()
        for drl in pairs.values():
            if len(drl) == 1:
                self.withoutDuplicates.append(drl[0])
            else:
                for d in drl:
#                    nTdebug("***%s" % d.format())
                    self.withDuplicates.append(d)
                #end for
            #end if
        #end for
#        nTdebug("For %s derived:\n%s" % (self, format(self)))
    #end def

    def format(self, allowHtml = False, showAll = False):  # pylint: disable=W0221
        a = len(self)
        b = len(self.withoutDuplicates)
        c = self.uniqueDistancesCount
        d = a - c
        e = c - b
        f = len(self.withDuplicates)
        if b + f != a:
            nTerror("Failed check in %s of b + f != a.")
        # end if
        msg = sprintf(
'''
classes
  intra-residual:     %4d
  sequential:         %4d
  medium-range:       %4d
  long-range:         %4d
  ambiguous:          %4d

counts
  total all:          %4d (A)        All DRs
  singly defined      %4d (B)        DRs for which no other DR has the same set of atom pairs.
  multiple atom pairs:%4d (F=A-B)
  multiple defined    %4d (E=C-B)
  total unique:       %4d (C=A-D)    Count of the set of unique atom pairs.
  duplicates:         %4d (D=A-C)

rmsd:              %7s +-%6s

violations
  <-0.1 A:            %4d (lower-bound violations, max:%5s)
  > 0.1 A:            %4d
  > 0.3 A:            %4d
  > 0.5 A:            %4d (upper-bound violations, max:%5s)

ROG score:         %7s
''',
                        len(self.intraResidual),
                        len(self.sequential),
                        len(self.mediumRange),
                        len(self.longRange),
                        len(self.ambiguous),

                        a, b, f, e, c, d,

                        val2Str(self.rmsdAv, "%7.3f", 7), val2Str(self.rmsdSd, "%6.3f", 6),
                        self.violCountLower, val2Str(self.violLowerMax, "%5.2f", 5),
                        self.violCount1, self.violCount3, self.violCount5, val2Str(self.violUpperMax, "%5.2f", 5),
                        self.rogScore
                      )
        if allowHtml:
            msg = addPreTagLines(msg)
        header = '%s DistanceRestraintList "%s" (%s,%d) %s\n' % (
            dots, self.name, self.status, len(self), dots)
        msg = header + msg
        msg += RestraintList.format(self, showAll = False)
        return msg
    #end def

#end class


class DihedralRestraint(Restraint):
    """
        DihedralRestraint class:

       GV 2 Oct 2007: removed residue and angle attributes.
       If the 4 atoms constitute a known dihedral angle, this can
       be retrieved with the retrieveDefinition method

       GV&AWSS: 10 Oct 2007, upper-limit adjustment
    """

#    project.valSets.DR_MAXALL_POOR = 3. # Normally 0.3 but set low for testing 1brv to
#    DR_MAXALL_BAD  = 5. # Normally 0.5 but set low for testing 1brv to
#    DR_THRESHOLD_OVER_BAD  = 0.3 # degrees.
#    DR_THRESHOLD_FRAC_BAD  = 0.5
#    DR_RMSALL_BAD  = 0.3 # Angstrom rms violations. # Normally 0.3 but set low for testing 1brv to

    def __init__(self, atoms, lower, upper, **kwds):
        if upper < lower:
            upper += 360.0
        Restraint.__init__(self,
                              lower = lower,
                              upper = upper,
                              **kwds
                       )
        self.dihedrals = NTlist() # list with dihedral values for each model
        self.cav = None      # Average dihedral value
        self.cv = None      # cv on dihedral

        self.setdefault('discontinuous', False)
        self.__CLASS__ = AC_LEVEL
        self.atoms = NTlist(*atoms)
        if None in self.atoms:
            self.isValid = False
    #end def

    def isValidForAquaExport(self):
        """Determine if the restraint can be exported to Aqua.
        Simplified to checking if there are 4 real atoms.
        """
        if not self.atoms:
#            nTdebug("Failed to find any atom in %s" % self)
            return False
        n = len(self.atoms)
        if n != 4:
#            nTdebug("Expected four atoms but found %d in:\n%s" % (n,self))
            return False
        for _i, atom in enumerate(self.atoms):
            if not atom:
#                nTdebug("Failed to find valid atom in:\n%s" % (i,self))
                return False
        return True
    #end def

    def criticize(self, project):
        """Only the self violations,violMax and violSd needs to be set before calling this routine"""
#        nTdebug( '%s (dih)' % self )
        self.rogScore.reset()

        if (project.valSets.AC_MAXALL_BAD != None) and (self.violMax >= project.valSets.AC_MAXALL_BAD):
            comment = 'violMax: %7.2f' % self.violMax
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_RED, comment)
        elif (project.valSets.AC_MAXALL_POOR != None) and (self.violMax >= project.valSets.AC_MAXALL_POOR):
            comment = 'violMax: %7.2f' % self.violMax
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_ORANGE, comment)

        if (project.valSets.AC_THRESHOLD_OVER_POOR != None) and (project.valSets.AC_THRESHOLD_FRAC_POOR != None):
            fractionAbove = getFractionAbove(self.violations, project.valSets.AC_THRESHOLD_OVER_POOR)
            if fractionAbove >= project.valSets.AC_THRESHOLD_FRAC_POOR:
                comment = 'fractionAbove: %7.2f' % fractionAbove
    #            nTdebug(comment)
                self.rogScore.setMaxColor(COLOR_ORANGE, comment)

        if (project.valSets.AC_THRESHOLD_OVER_BAD != None) and (project.valSets.AC_THRESHOLD_FRAC_BAD != None):
            fractionAbove = getFractionAbove(self.violations, project.valSets.AC_THRESHOLD_OVER_BAD)
            if fractionAbove >= project.valSets.AC_THRESHOLD_FRAC_BAD:
                comment = 'fractionAbove: %7.2f' % fractionAbove
    #            nTdebug(comment)
                self.rogScore.setMaxColor(COLOR_RED, comment)

        if (project.valSets.AC_RMSALL_BAD != None) and (self.violSd >= project.valSets.AC_RMSALL_BAD):
            comment = 'violSd: %7.2f' % self.violSd
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_RED, comment)
        elif (project.valSets.AC_RMSALL_POOR != None) and (self.violSd >= project.valSets.AC_RMSALL_POOR):
            comment = 'violSd: %7.2f' % self.violSd
#            nTdebug(comment)
            self.rogScore.setMaxColor(COLOR_ORANGE, comment)

        if project.valSets.FLAG_MISSING_COOR:
            #modelCount = self.getModelCount()
                #atms = atm.realAtoms()
                #for a in atms:
            for a in self.atoms:                # GV says: no realAtoms should be called; these are dihedrals
                                                # and should be properly defined by their \atoms.

                #if len(a.coordinates) < modelCount:
                if a and a.hasMissingCoordinates(): # gv has mase this into a method because the getModelCount()
                    # can crash when reading back NRG dataset because of their
                    # incompleteness
                    msg = "Missing coordinates in dihedral (%s)" % a.toString()
#                    nTdebug(msg)
                    self.rogScore.setMaxColor(COLOR_RED, msg)
    #end def

    def isChi2TyrOrPhe(self):
        lastAtom = self.atoms[3]
        atomName = getDeepByKeysOrAttributes(lastAtom,'name')
        if atomName != 'CD1':
            return False
        lastAtomRes = lastAtom._parent
        lastAtomResType = getDeepByKeysOrAttributes(lastAtomRes, 'db', 'name')
        if lastAtomResType in ['TYR', 'PHE']:
            return True
        return False

    def calculateAverage(self):
        """Calculate the values and violations for each model
        return cav and cv tuple or (None, None) tuple on error
        """
        errorExit = (None, None)
        if len(self.atoms) != 4 or (None in self.atoms):
            nTerror('DihedralRestraint: invalid dihedral definition %s', self.atoms)
            return errorExit
        #end if

        if None in self.atoms.zap('meanCoordinate'):
            nTerror('DihedralRestraint: atom(s) without coordinates %s', self.atoms)
            return errorExit
        #end if

#        coorList = self.atoms.zap('coordinates')
#        if len( coorList ) == 0:
#            nTerror('DihedralRestraint: atom(s) without any coordinates %s', self.atoms)
#            return (None, None)
#        #end if

        modelCount = self.getModelCount()
        if modelCount == 0:
            nTerror('DihedralRestraint: no structure models')
            return errorExit
        #end if
#        lenCoorListExpected = 4 * modelCount
#        if len( coorList ) != lenCoorListExpected:
#            nTerror('DihedralRestraint: atom(s) without all coordinates %s', self.atoms)
#            return (None, None)
#        #end if


        #set the default values (JFD: this needs to be fully done in initializer in case code fails as for issue 222)
        self.dihedrals = NTlist() # list with dihedral values for each model
        self.cav = None      # Average dihedral value
        self.cv = None      # cv on dihedral

        self.violations = NTlist() # list with violations for each model
        self.violCount1 = 0        # Number of violations over 1 degree
        self.violCount3 = 0        # Number of violations over 3 degrees
        self.violCount5 = 0        # Number of violations over 5 degrees
        self.violMax = 0.0      # Maximum violation
        self.violAv = 0.0      # Average violation
        self.violSd = 0.0      # Sd of violations

        #find the range to store these dihedral values
        plotpars = plotParameters.getdefault(self.retrieveDefinition()[1], 'dihedralDefault')
        considerSymmetry = self.isChi2TyrOrPhe() # Hack for Phe/Tyr CHI2
        lastAtom = self.atoms[3]
        ssaPartner = None
        if considerSymmetry:
#            ssaPartner = lastAtom.getStereoPartner()
            try:
                ssaPartner = lastAtom._parent.CD2
            except:
                pass
#            nTdebug("ssaPartner: %s" % ssaPartner)
            if ssaPartner != None:
                considerSymmetry = True
            else:
                nTwarning("DihedralRestraint: no lastAtom's ssa for %s so ignoring symmetry on violation." % self)
                considerSymmetry = False

        if considerSymmetry:
            jLoopList = [ lastAtom, ssaPartner ]
        else:
            jLoopList = [ lastAtom ]

        try:
            # For each model we'll use the atom HD1 or HD2 that has the smallest violation or HD1 if neither one
            # is violated.
            for i in range(modelCount):
                dList = []
                vList = []
                for _j1, lastAtom2 in enumerate(jLoopList):
#                    nTdebug('i, _j1, lastAtom2, considerSymmetry: %s %s %s %s' % (i,_j1,lastAtom2, considerSymmetry))
                    atomList = [self.atoms[k] for k in range(3)]
                    atomList.append( lastAtom2 )
                    coorList = [ atom.coordinates[i] for atom in atomList]
                    d = nTdihedralOpt( *coorList )
                    if d == None:
#                        nTdebug("Failed to calculate an angle; which can happen if a coordinate is missing.")
                        continue
                    dList.append( d )
                # end for _j1
                nTlimit(dList, plotpars.min, plotpars.max)
                for _j2 in range(len(dList)):
                    v = violationAngle(value = dList[_j2], lowerBound = self.lower, upperBound = self.upper)
                    if v == None:
                        nTwarning("Failed to calculate a violation angle.")
                        return errorExit
                    vList.append( v )
                # end for _j2
                jSelected = 0
                if considerSymmetry:
                    fvList = [ math.fabs(x) for x in vList]
                    if len(fvList) == 2:
                        if fvList[1] < fvList[0]:
                            jSelected = 1
#                    nTdebug("Comparing fviolations for %s %s" % ( self, fvList))
                # end if
#                nTdebug("Comparing distances for %s %s" % ( self, dList))
#                nTdebug("Comparing violations for %s %s" % ( self, vList))
#                nTdebug("jSelected %s" % jSelected)
                self.dihedrals.append(dList[jSelected])
                self.violations.append(vList[jSelected])
#                nTdebug("self.dihedrals %s" % self.dihedrals)
#                nTdebug("self.violations %s" % self.violations)

                fv = math.fabs(vList[jSelected])
                if fv > 1.0:
                    self.violCount1 += 1
                if fv > 3.0:
                    self.violCount3 += 1
                if fv > 5.0:
                    self.violCount5 += 1
                if fv > self.violMax:
                    self.violMax = fv
                #end if
            #end for all models
        except:
#            NTtracebackError() # DEFAULT this is disabled.
#            nTdebug("Ignoring violations for %s" % self.format() )
            pass # ignore missing coordinates. They're reported by criticize()

        self.violAv, self.violSd, _n = self.violations.average()
        # The CV is hard to calculate for the symmetry case detailed above. TODO:
        self.cav, self.cv, _n = self.dihedrals.cAverage(plotpars.min, plotpars.max)
        return(self.cav, self.cv)
    #end def

    def retrieveDefinition(self):
        """
        Retrieve a (<Residue>, angleName, <AngleDef>) tuple from
        the molecule._dihedralDict
        or
        (None,None,None) on error
        """
        resultError = (None, None, None)
        if (not self.atoms or (None in self.atoms)):
            return resultError
        #end if

#        mol = self.atoms[0].residue.chain.molecule
        mol = self.atoms[0].getMolecule()
        if mol == None:
            return resultError
        if mol._dihedralDict.has_key(tuple(self.atoms)):
            return mol._dihedralDict[tuple(self.atoms)]
        return resultError
        #end if
    #end def

    def getName(self):
        """
        Construct a name;
        have to do dynamically because upon restoring, the atoms are not yet defined
        """
        res, name, _tmp = self.retrieveDefinition()
        if res and name:
            return res.name + '.' + name
        else:
            return ''
    #end def

    def getDihedralName(self):
        """
        Construct a name;
        have to do dynamically because upon restoring, the atoms are not yet defined
        """
        _res, name, _tmp = self.retrieveDefinition()
        return name
    #end def

    def format(self):  # pylint: disable=W0221
        # set the last string to something readable in terms of known dihedrals, or just the atoms if nothing is found
        s = self.getName()
        if len(s) == 0:
            s = self.atoms.format('%-11s ')
        return  \
            sprintf('%-25s %-6s (Target: %s %s)  (Models: cav %6s cv %7s)  ' + \
                    '(Violations: av %4s max %4.1f counts %2d,%2d,%2d) %s',
                     self, self.rogScore,
                     val2Str(self.lower, "%6.1f", 6), # GV: does not fit in 4 fields, i.e -160.1
                     val2Str(self.upper, "%6.1f", 6),
                     val2Str(self.cav, "%6.1f", 6),
                     val2Str(self.cv, "%7.4f", 7),
                     val2Str(self.violAv, "%4.1f", 4),
                     self.violMax, self.violCount1, self.violCount3, self.violCount5,
                     s
                    )
    #end def
#end class

# Too many ancestors (8/7) pylint: disable=R0901
class DihedralRestraintList(RestraintList):
    'List of dihedral angle restraints.'
#    export2cyana = exportDihedralRestraint2cyana

    def __init__(self, name, status = 'keep'):
        RestraintList.__init__(self, name, status = status)
        self.__CLASS__ = ACL_LEVEL
    #end def

    def criticize(self, project, toFile = True):
        """
        Criticize restraints of this list; infer own ROG score from individual restraints.
        """
#        nTdebug('DihedralRestraintList.criticize %s', self)
        self.rogScore.reset()

        for dr in self:
            dr.criticize(project)
            self.rogScore.setMaxColor(dr.rogScore.colorLabel, comment = 'Cascaded from: %s' % dr)
        if toFile:
            path = project.moleculePath('analysis', self.name + '.txt')
            f = file(path, 'w')
            fprintf(f, '%s\n\n', self.format())
            for dr in self:
                fprintf(f, '%s\n', dr.format())
            #end for
            f.close()
            nTdetail('==> Analyzing %s, output to %s', self, path)
        #end if
    #end def

    def analyze(self, calculateFirst = True):
        """
        Calculate averages for every restraint.
        Return <rmsd>, sd and total violations over 1, 3, and 5 degrees as tuple
        or (None, None, None, None, None) on error
        """
        errorResult = (None, None, None, None, None)
#        nTdebug('DihedralRestraintList.analyze: %s', self)

        if not len(self):
            nTwarning('DihedralRestraintList.analyze: "%s" empty list', self.name)
            return errorResult
        #end if

        modelCount = self.getModelCount()
        if not modelCount:
            nTerror('DihedralRestraintList.analyze: "%s" modelCount 0', self.name)
            return errorResult
        #end if

        self.rmsd = nTfill(0.0, modelCount)
        self.violCount1 = 0
        self.violCount3 = 0
        self.violCount5 = 0
        for dr in self:
            if calculateFirst:
                (cav, _cv) = dr.calculateAverage()
                if cav == None:
#                    nTdebug("Failed to calculate average for: " + self.format())
                    continue # skipping dihedral with a missing coordinate or so.
            self.violCount1 += dr.violCount1
            self.violCount3 += dr.violCount3
            self.violCount5 += dr.violCount5

            countDrViolations = len(dr.violations)
            if countDrViolations > modelCount:
                nTcodeerror("Found more violations (%s) for this restraint (%s) than models (%s)" % ( countDrViolations, dr, modelCount))
                return errorResult

            for i in range(countDrViolations): # happened in entry 1bn0 that violations were not defined.
#            for i in range(0, modelCount):
                v = dr.violations[i]
                self.rmsd[i] += v * v
            #end for
        #end for

        for i in range(0, modelCount):
            self.rmsd[i] = math.sqrt(self.rmsd[i] / len(self))
        #end for

        self.rmsdAv, self.rmsdSd, _n = nTaverage(self.rmsd)
        return (self.rmsdAv, self.rmsdSd, self.violCount1, self.violCount3, self.violCount5)
    #end def

    def format(self, allowHtml = False):
        msg = sprintf(
'''
rmsd:                 %7s +-%6s
violations >1 degree: %4d
violations >3 degree: %4d
violations >5 degree: %4d

ROG score:            %s
''', val2Str(self.rmsdAv, "%7.3f", 7),
                        val2Str(self.rmsdSd, "%6.3f", 6),
                        self.violCount1,
                        self.violCount3,
                        self.violCount5,
                        self.rogScore
                      )
        if allowHtml:
            msg = addPreTagLines(msg)
        header = '%s DihedralRestraintList "%s" (%s,%d) %s\n' % (
            dots, self.name, self.status, len(self), dots)
        msg = header + msg
        return msg
    #end def

    def formatHtml(self, allowHtml = False):
        header = self.name
        if hasattr(self, 'rogScore'):
            if self.rogScore.isCritiqued():
                header = '<font color="%s">%s</font>' % (self.rogScore.colorLabel, header)
        header = '<h3>DihedralRestraintList %s</h3>' % header

        msg = '''%s
<BR>
<table>
<TR><TD>rmsd:               </TD><TD> %s +- %s                    </TD></TR>
<TR><TD>violations > 1 degree</TD><TD align="right"> %4d                          </TD></TR>
<TR><TD>violations > 3 degrees</TD><TD align="right"> %4d                          </TD></TR>
<TR><TD>violations > 5 degrees</TD><TD align="right"> %4d                          </TD></TR>
</table>
''' % (
    header,
    val2Str(self.rmsdAv, "%7.3f", 7), val2Str(self.rmsdSd, "%6.3f", 6),
    self.violCount1,
    self.violCount3,
    self.violCount5
  )
        return msg
    #end def

    def isFromTalos(self):
        return self.name.startswith(TALOSPLUS_LIST_STR)
#end class


class RDCRestraint(DistanceRestraint):
    """Residual Dipolar Coupling restraint
    much like DistanceRestraint having also atomPairs.
    """
    def __init__(self, atomPairs = [], lower = 0.0, upper = 0.0, **kwds):

        DistanceRestraint.__init__(self,
                              atomPairs = atomPairs,
                              lower = lower,
                              upper = upper,
                              **kwds
                       )
        self.__CLASS__ = RDC_LEVEL
        self.rdcs = None     # list with backcalculated rdc values for each model, None indicates no analysis done
        # copied from dihedral; to be implemented
        self.cav = None     # Average dihedral value
        self.cv = None     # cv on dihedral
    #end def


    def calculateAverage(self):
        """Calculate the values and violations for each model
           return cav and cv tuple or (None, None) tuple on error
        """

        modelCount = self.getModelCount()
        if not modelCount:
            nTerror('Error RDCRestraint: no structure models\n')
            return (None, None)
        #end if

#TODO needs work
#        if len(self.atoms) != 2 or None in self.atoms:
#            nTerror('Error RDCRestraint: invalid rdc definition %s\n', self.atoms )
#            return (None,None)
#        #end if
#
#        if None in self.atoms.zap('meanCoordinate'):
#            nTerror('Error RDCRestraint: atom without coordinates %s\n', self.atoms )
#            return (None,None)
#        #end if

        #set the default values; do not remove or it will crash upon restoring RDC lists
        self.rdcs = nTfill(0.0, modelCount) # list with dihedral values for each model
        self.cav = NaN      # Average dihedral value
        self.cv = NaN      # cv on dihedral
#        self.min        = 0.0      # Minimum dihedral
#        self.max        = 0.0      # Max dihedral
        self.violations = nTfill(0.0, modelCount) # list with violations for each model
        self.violCount1 = 0        # Number of violations over 1 degree
        self.violCount3 = 0        # Number of violations over 3 degrees
        self.violCount5 = 0        # Number of violations over 5 degrees
        self.violMax = 0.0      # Maximum violation
        self.violAv = 0.0      # Average violation
        self.violSd = 0.0      # Sd of violations

        return(None, None)
    #end def

    def _names(self):
        """
        Internal routine: generate string from atomPairs
        """
        s = ''
        for p in self.atomPairs:
            s = s + sprintf('(%-11s - %11s)   ', p[0].cName(1), p[1].cName(1))
        #end for
        return s.strip()
    #end def

    def format(self): # pylint: disable=W0221
#        s = '('
#        for p in self.atoms:
#            s = s + sprintf('%-11s ', p.cName(1) )
#        #end for
#        s = s.strip() + ')'
        return  sprintf('%-25s %-6s (Target: %6.1f %6.1f) %s',
                        str(self), self.rogScore,
                        self.lower, self.upper,
                        self._names()
                       )

    #end def
#end class

# Too many ancestors (8/7) pylint: disable=R0901
class RDCRestraintList(RestraintList):
    """List of RDCRestraint"""

    def __init__(self, name, status = 'keep'):
        RestraintList.__init__(self, name, status = status)
        self.__CLASS__ = RDCL_LEVEL
    #end def

    def analyze(self, calculateFirst = True):
        """
        Calculate averages for every restraint.

        """

#        nTdebug('RDCRestraintList.analyze: %s', self)

        if len(self) == 0:
            nTerror('RDCRestraintList.analyze: "%s" empty list', self.name)
            return (None, None, None, None, None)
        #end if

        modelCount = 0
        firstRestraint = self[0]
        if not hasattr(firstRestraint, "atoms"):
            nTwarning("Failed to get the model count for no atoms are available in the first RDC restraint.")
            nTwarning("See also issue: %s%d" % (issueListUrl, 133))
        else:
            if len(self[0].atomPairs):
                modelCount = self[0].atomPairs[0][0].residue.chain.molecule.modelCount
        #end if

        if not modelCount: # JFD notes eg reading $CINGROOT/Tests/data/ccpn/2hgh.tgz
            nTerror('RDCRestraintList.analyze: "%s" modelCount 0', self.name)
            return (None, None, None, None, None)
        #end if

        self.rmsd = nTfill(0.0, modelCount)
        self.violCount1 = 0
        self.violCount3 = 0
        self.violCount5 = 0
        for dr in self:
            if calculateFirst:
                dr.calculateAverage()
            if dr.rdcs == None or dr.violations == None:
                nTerror('RDCRestraintList.analyze: skipping restraint %s', dr)
            else:
                self.violCount1 += dr.violCount1
                self.violCount3 += dr.violCount3
                self.violCount5 += dr.violCount5
                for i in range(0, modelCount):
                    self.rmsd[i] += dr.violations[i] * dr.violations[i]
                #end for
            #end if
        #end for

        for i in range(0, modelCount):
            self.rmsd[i] = math.sqrt(self.rmsd[i] / len(self))
        #end for

        self.rmsdAv, self.rmsdSd, dummy_n = nTaverage(self.rmsd)
        return (self.rmsdAv, self.rmsdSd, self.violCount1, self.violCount3, self.violCount5)
    #end def

    def criticize(self, project, toFile = True):
        """
        Criticize restraints of this list; infer own ROG score from individual restraints.

        TODO: Need implementation
        """
#        nTdebug('RDCRestraintList.criticize %s', self)
        self.rogScore.reset()

#        for dr in self:
#            dr.criticize(project)
#            self.rogScore.setMaxColor( dr.rogScore.colorLabel, comment='Cascaded from: %s' % dr )
#        if toFile:
#            path = project.moleculePath('analysis', self.name + '.txt')
#            f = file( path, 'w')
#            fprintf(f, '%s\n\n', self.format())
#            for dr in self:
#                fprintf(f, '%s\n', dr.format())
#            #end for
#            f.close()
#            nTdetail('Distance restraint analysis %s, output to %s', self, path)
#        #end if
    #end def


    def format(self, allowHtml = False):
        msg = sprintf('rmsd: %7.3f %6.3f',
                      self.rmsdAv, self.rmsdSd)
        if allowHtml:
            msg = addPreTagLines(msg)
        header = '%s RDCRestraintList "%s" (%s,%d) %s\n' % (
            dots, self.name, self.status, len(self), dots)
        msg = header + msg
        return msg
    #end def

#    def formatHtml( self ):
##TODO clean up
#        header = self.name
#        if hasattr(self, 'rogScore'):
#            if self.rogScore.isCritiqued():
#                header = '<font color="%s">%s</font>' % (self.rogScore.colorLabel, header)
#        header = '<h3>RDCRestraintList %s</h3>' % header
#
#        msg = '''%s
#<BR>
#<table>
#<TR><TD>rmsd</TD>               <TD> %s +- %s                    </TD></TR>
#<TR><TD>violations > 1 </TD><TD align="right"> %4d                          </TD></TR>
#<TR><TD>violations > 3 </TD><TD align="right"> %4d                          </TD></TR>
#<TR><TD>violations > 5 </TD><TD align="right"> %4d                          </TD></TR>
#</table>
#''' % (
#    header,
#    val2Str(self.rmsdAv, "%7.3f", 7), val2Str(self.rmsdSd, "%6.3f", 6),
#    self.violCount1,
#    self.violCount3,
#    self.violCount5
#  )
#        return msg
#    #end def

    def formatHtml(self, allowHtml = False):
        header = self.name
        if hasattr(self, 'rogScore'):
            if self.rogScore.isCritiqued():
                header = '<font color="%s">%s</font>' % (self.rogScore.colorLabel, header)
        header = '<h3>RDCRestraintList %s</h3>' % header

        msg = '''%s
<BR>
<table>
<TR><TD>count               </TD><TD align="right">%4d</TD></TR>
</table>
''' % (
    header,
    len(self)
  )
        return msg
    #end def
#end class






def path(*args):
    """
    Return a path from arguments relative to cing root
    """
    return os.path.join(cingRoot, *args)
#end def

def shift(atm):
    return atm.shift()
#end def


def getFractionAbove(valueList, threshold):
    """Return the fraction that is above the threshold. Return .0 if
    list contains no values.
    """
    if not valueList:
        return .0
    n = 0.0 # enforce float arithmics for below division.
    for v in valueList:
        if not v: # catch None and pure zero
            continue
        if v > threshold:
            n += 1.
    return n / len(valueList)

