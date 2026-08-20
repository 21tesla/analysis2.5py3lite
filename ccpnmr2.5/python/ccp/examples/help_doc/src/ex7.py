import os
import Tkinter

from memops.api.Implementation import MemopsRoot

from ccpnmr.format.converters.CnsFormat import CnsFormat

if __name__ == '__main__':

    # CCPN project.
    project = MemopsRoot(name = 'readCns')

    guiRoot = Tkinter.Tk()

    # CNS object.
    cnsObj = CnsFormat(project, guiRoot)

    # Read in a sequence - this will create the molecular system with
    # all the atom information.

    # Note, many of the popups can be avoided when the right information
    # is passed in (see the readSequence() function in the
    # ccpnmr.format.converters.DataFormat class.

    # Set the location of the CNS files.
    cnsDir = '../data/cns'

    # Read in a sequence to set molecule and molSystem objects.
    cnsCoord     = 'cns_1.pdb'
    cnsCoordFile = os.path.join(cnsDir, cnsCoord)

    # readSequence() returns a list of chains.
    chains = cnsObj.readSequence(cnsCoordFile,
                                 minimalPrompts = 1)

    # Read in a distance constraint list.
    distConst     = 'n15noesy.tbl'
    distConstFile = os.path.join(cnsDir, distConst)

    constList = cnsObj.readDistanceConstraints(distConstFile,
                                               minimalPrompts = 1)

    # Do some preliminary data model navigation to get input parameters for
    # linkResonances.

    # An nmrConstraintStore links a group of constraint files.
    # A structureGeneration links an nmrConstraintStore with a set 
    # of structures.
    nmrConstStore = constList.nmrConstraintStore
    strucGen      = nmrConstStore.findFirstStructureGeneration()

    # Run linkResonances (this will generate a lot of output to the shell).

    # The current options are the 'safest' to maintain the original
    # information, although bear in mind that all atoms in the original
    # list are considered to be stereospecifically assigned.

    # Set forceDefaultChainMapping to 0 if you want to interactively link the
    # chains in the CCPN data model to the information from the constraint file.
    cnsObj.linkResonances(forceDefaultChainMapping = 1,
                          globalStereoAssign = 1,
                          setSingleProchiral = 1,
                          setSinglePossEquiv = 1,
                          strucGen = strucGen,
                          minimalPrompts = 0)

    # Navigate the data model to get a list of atoms per constraint item.
    for distConst in constList.sortedConstraints():

        print 'Constraint [%d]: [%.1f] - [%.1f]' % (
            distConst.serial, distConst.lowerLimit, distConst.upperLimit)

        for constItem in distConst.sortedItems():

            # Now list the atoms linked to each of the two resonances
            # associated with this item.
            atomList = []

            # Sometimes there may also be an ordered<Class> method.
            resonanceList = constItem.orderedResonances

            # Otherwise, use the usual sorted<Class> method.
            if not resonanceList:
                resonanceList = constItem.sortedResonances()

            # Here, resonanceList should always have 2 resonances.
            assert(len(resonanceList) == 2)

            for resonance in resonanceList:
                resAtomList = []

                if resonance.resonanceSet:
                    resonanceSet = resonance.resonanceSet

                    for atomSet in resonanceSet.sortedAtomSets():
                        for atom in atomSet.sortedAtoms():
                            resAtomList.append('%d.%s' % (
                                atom.residue.seqCode, atom.name) )

                resAtomList.sort()
                resAtomString = ','.join(resAtomList)

                atomList.append(resAtomString)

            print '  [%s] - [%s]' % (atomList[0], atomList[1])

        print
