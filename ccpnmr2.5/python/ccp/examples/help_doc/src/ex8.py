import os
import Tkinter

from memops.api.Implementation import MemopsRoot

from ccpnmr.format.converters.NmrViewFormat import NmrViewFormat
from ccpnmr.format.converters.XEasyFormat import XEasyFormat

# Pre-defined functions to create an NMR experiment and data source.
from ccpnmr.format.general.Util import (
    createExperiment, getRefExpFromOldExpType, createPpmFreqDataSource)

if __name__ == '__main__':

    # CCPN project.
    project = MemopsRoot(name = 'nmrView2XEasy')

    gui = Tkinter.Tk()

    # NmrView object.
    nmrViewObj = NmrViewFormat(project, gui)

    # Location of the NmrView files.
    nmrViewDir = '../data/nmrView'

    # Read in a sequence to set molecule and molSystem objects.
    nmrViewSeq     = 'nmrView.seq'
    nmrViewSeqFile = os.path.join(nmrViewDir, nmrViewSeq)

    nmrViewObj.readSequence(nmrViewSeqFile, minimalPrompts = 1)

    # At this point we need an nmrProject object to save chemical shifts,
    # peak lists, etc.
    nmrProject = project.currentNmrProject = \
        project.newNmrProject(name = project.name)

    # Read in a file containing chemical shifts.
    nmrViewChemShifts    = 'ppm.out'
    nmrViewChemShiftFile = os.path.join(nmrViewDir, nmrViewChemShifts)

    nmrViewObj.readShifts(nmrViewChemShiftFile, minimalPrompts = 1)

    # Find and print this shiftList in the CCPN data model.
    shiftList = nmrProject.findFirstMeasurementList(className = 'ShiftList')

    print 'ShiftList: [%s]' % shiftList

    # Create a reference NMR experiment. This is based on the 
    # NmrExpPrototype setup.
    refExpType = 'noesy_hsqc_HNH.hhn'  # 'Old' experiment type name.
    refExp = getRefExpFromOldExpType(project, refExpType)

    # Should print 'H_H[N].NOESY'.
    print 'New experiment type name: [%s]' % refExp.name

    # Create a real NMR experiment from this object and link it to the
    # CCPN project.
    keywds = {'expName': 'noesyTest', 'refExperiment': refExp}
    nmrExp = createExperiment(project, **keywds)

    # Create a dataSource object (with the right number of dimensions etc.)
    # to hold peak list information.
    keywds2 = {'dsName': 'test', 'dsType': 'processed', 'numDim': nmrExp.numDim}
    nmrDataSource = createPpmFreqDataSource(nmrExp, **keywds2)

    # Read in a peak list connected to the experiment/dataSource that 
    # was just created.
    nmrViewPeakLists    = 'nmrView.xpk'
    nmrViewPeakListFile = os.path.join(nmrViewDir, nmrViewPeakLists)

    # Do the reading - this will throw up a Tkinter popup window.
    nmrViewObj.readPeaks(nmrViewPeakListFile,
                         dataSource = nmrDataSource,
                         minimalPrompts = 1)

    # Find and print this peakList in the CCPN data model.
    peakLists = nmrDataSource.sortedPeakLists()
    print 'Peaklists: [%s]' % peakLists

    # Run linkResonances.

    # This is an important sub-routine that links together the NMR data
    # (peaks, shifts), which are found in the 'Resonance' objects, to the
    # sequence information (molecules, chains, residues, atoms). Without doing
    # this, these objects would not be linked to each other in the CCPN project.

    # In basic terms, linkResonances links up the NMR information to the
    # atom information - hence, assigning the NMR information.

    # Many options are available - see ccpnmr/format/process/linkResonances.py.
    # Also see what happens when you comment this function out.
    nmrViewObj.linkResonances(setSingleProchiral = 0,
                              setSinglePossEquiv = 0,
                              minimalPrompts = 1)

    # Note that at this stage everything is inside the data model. Any code that
    # works with the data model can be run at this stage (e.g., creating 
    # a chemical shift list from a set of peak assignments).

    # Location to save XEasy files.
    curDir   = os.path.abspath('../data')
    xeasyDir = os.path.join(curDir, 'xEasy')

    if not os.path.exists(xeasyDir):
        os.mkdir(xeasyDir)

    # XEasy object.
    xeasyObj = XEasyFormat(project, gui)

    # Write out sequence in XEasy format.
    chains = project.sortedMolSystems()[0].sortedChains()
 
    xeasySeq     = 'xeasy.seq'
    xeasySeqFile = os.path.join(xeasyDir, xeasySeq)

    # Use the chains object to write out the sequence to a file.
    xeasyObj.writeSequence(xeasySeqFile, chains = chains, minimalPrompts = 1)

    # Write out chemical shifts in XEasy format.
    xeasyChemShifts    = 'xeasy.prot'
    xeasyChemShiftFile = os.path.join(xeasyDir, xeasyChemShifts)

    # Use the shiftList object made earlier for printing.
    xeasyObj.writeShifts(xeasyChemShiftFile,
                         measurementList = shiftList,
                         minimalPrompts = 1)

    # Write out peak lists in XEasy format.
    xeasyPeakLists    = 'xeasy.xpk'
    xeasyPeakListFile = os.path.join(xeasyDir, xeasyPeakLists)

    # Use the peakLists object made earlier for printing - throws up a popup.
    xeasyObj.writePeaks(xeasyPeakListFile,
                        peakLists = peakLists,
                        minimalPrompts = 1)

    # Save the project.
    project.saveModified()

    # Some data model navigation examples.
    allResonances = project.currentNmrProject.sortedResonances()
    numResonances = len(allResonances)

    print '\nNumber of resonances: [%d]' % numResonances

    # Information about the NMR experiments in this project.
    allExpts = project.currentNmrProject.sortedExperiments()

    for nmrExp in allExpts:
        print 'Experiment: [%d], called [%s], type [%s], with [%d] dims.' % (
            nmrExp.serial, nmrExp.name, nmrExp.experimentType, nmrExp.numDim)

        # This line prints a different name as this is a 'reverse' experiment.
        print 'Reference pathway: [%s]' % (
            nmrExp.refExperiment.nmrExpPrototype.name)

    # List to save all the residues that are associated with these resonances.
    seqCodes = []

    for resonance in allResonances:
        if resonance.resonanceSet:
            resonanceSet = resonance.resonanceSet

            for atomSet in resonanceSet.sortedAtomSets():
                seqCode = atomSet.sortedAtoms()[0].residue.seqCode

                if seqCode not in seqCodes:
                    seqCodes.append(seqCode)

    seqCodes.sort()

    print 'Residue sequence codes connected to resonances:'
    print seqCodes, '\n'
