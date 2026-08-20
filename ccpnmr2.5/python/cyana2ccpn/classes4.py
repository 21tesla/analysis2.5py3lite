from cing.Libs.NTutils import *
from cing.core.classes2 import RestraintList

class Restraint(NTdict):
    """
    Super class for DistanceRestraint etc.
    On initialization the atom pairs will be tested for validity and reported in self.isValid
    """
    def __init__(self, value, error, weight, **kwds):
        NTdict.__init__(self,value=value,
                              error=error,
                              weight=weight,
                              **kwds
                          )
        self.__CLASS__ = None     # Set by sub class.

    #end def

class RDCRestraint(Restraint):
    """Residual Dipolar Coupling restraint
    much like DistanceRestraint having also atomPairs.
    """
    def __init__(self, atomPairs = NTlist(), value = 0.0, error = 0.0, weight = 0.0, **kwds):



        Restraint.__init__(self,
                              atomPairs = atomPairs,
                              value=value,
                              error=error,
                              weight=weight,
                              **kwds
                       )
        self.__CLASS__ = RDC_LEVEL

        self.atomPairs = NTlist()
        print self.atomPairs()

    def appendPair(self, pair):
        """ pair is a (atom1,atom2) tuple

        check if atom1 already present, keep order
        otherwise: keep atom with lower residue index first
        Return True on error.
        """
        # GV says; order needs to stay: is being used for easier
        # (manual) analysis.


        if pair[0] == None or pair[1] == None:
            nTerror('RDCRestraint.appendPair: invalid pair %s', str(pair))
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

    #end def

class RDCRestraintList(RestraintList):
    """
    Class based on NTlist that holds distanceRestraints.
    Also manages the "id's".
    Sort  by item of DistanceRestraint Class
    """
    # use the same spelling through out.
    def __init__(self, name, status = 'keep'):
        RestraintList.__init__(self, name, status = status)

class ShiftRestraint( NTdict ):
    """Residual Dipolar Coupling restraint
    much like DistanceRestraint having also atomPairs.
    """

    def __init__(self, shift, error, **kwds):
        NTdict.__init__(self, shift = shift,
                              error = error,
                              **kwds
                       )
        self.__CLASS__ = None     # Set by sub class.

    #end def

class ChemicalShiftRestraint(ShiftRestraint):
    """ChemicalShiftRestraint class:
       atom: atom,
       shift and error bounds.
    """
    def __init__(self, atom = None, shift = None, error = None, **kwds):
        ShiftRestraint.__init__(self, shift = shift, error = error, **kwds)
        self.__CLASS__ = DR_LEVEL
        self.atom = atom

class ChemicalShiftRestraintList(RestraintList):
    """
    Class based on NTlist that holds distanceRestraints.
    Also manages the "id's".
    Sort  by item of DistanceRestraint Class
    """
    # use the same spelling through out.
    def __init__(self, name, status = 'keep'):
        RestraintList.__init__(self, name, status = status)
