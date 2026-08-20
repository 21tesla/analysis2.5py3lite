import os
import Tkinter

from memops.api.Implementation import MemopsRoot

from ccpnmr.format.converters.NmrViewFormat import NmrViewFormat

from ccpnmr.format.general.Util import (
    createExperiment, getRefExpFromOldExpType, createPpmFreqDataSource)

# Specific class not in Format Converter used to export NMR-STAR files.
from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport

if __name__ == '__main__':

    # This code follows on from ex8.py, but writes the project
    # to NMR-STAR format.
    project = MemopsRoot(name = 'nmrStarTest')
    gui = Tkinter.Tk()

    nmrViewObj = NmrViewFormat(project, gui)

    nmrViewDir = '../data/nmrView'

    nmrViewSeq     = 'nmrView.seq'
    nmrViewSeqFile = os.path.join(nmrViewDir, nmrViewSeq)

    nmrViewObj.readSequence(nmrViewSeqFile, minimalPrompts = 1)

    nmrProject = project.currentNmrProject = \
        project.newNmrProject(name = project.name)

    nmrViewChemShifts    = 'ppm.out'
    nmrViewChemShiftFile = os.path.join(nmrViewDir, nmrViewChemShifts)

    nmrViewObj.readShifts(nmrViewChemShiftFile, minimalPrompts = 1)

    refExpType = 'noesy_hsqc_HNH.hhn'
    refExp = getRefExpFromOldExpType(project, refExpType)
    nmrExp = createExperiment(project, 'noesyTest', refExp)

    # Find the molSystem object in the CCPN project.
    molSystem = project.sortedMolSystems()[0]

    # Make a link between the molSystem and the experiment objects.
    nmrExp.addMolSystem(molSystem)

    nmrDataSource = createPpmFreqDataSource(nmrExp, 'test', 'processed',
                                            nmrExp.numDim)

    nmrViewPeakLists    = 'nmrView.xpk'
    nmrViewPeakListFile = os.path.join(nmrViewDir, nmrViewPeakLists)

    nmrViewObj.readPeaks(nmrViewPeakListFile,
                         dataSource = nmrDataSource,
                         minimalPrompts = 1)

    nmrViewObj.linkResonances(setSingleProchiral = 0,
                              setSinglePossEquiv = 0,
                              minimalPrompts = 1)

    # New part of the code - need a BMRB Entry object for writing to
    # an NMR-STAR file.
    entryStore = project.newNmrEntryStore(name = project.name)
    entry      = entryStore.newEntry(name = project.name)

    # Set links from the objects created earlier to this Entry object.
    entry.setMolSystem(molSystem)
    entry.addExperiment(nmrExp)

    # Get the shiftList object containing the chemical shifts.
    shiftList = nmrProject.findFirstMeasurementList(className = 'ShiftList')
    entry.addMeasurementList(shiftList)

    # Peak lists are not yet in the NmrStarExport dictionary.
    # This dictionary can be found in pdbe/nmrStar/IO/Ccpn_To_NmrStar.py.
    peakList = \
  entry.sortedExperiments()[0].sortedDataSources()[0].sortedPeakLists()[0]

    print 'Peak list: [%s]' % peakList

    # Directory and file to save the NMR-STAR file in.
    curDir     = os.path.abspath('../data')
    nmrStarDir = os.path.join(curDir, 'nmrStar')

    if not os.path.exists(nmrStarDir):
        os.mkdir(nmrStarDir)

    # Name of the NMR-STAR file.
    outNmrStarFile = os.path.join(nmrStarDir, 'nmrStar.str')

    # Code to create an NMR-STAR version 3.1 file from a BMRB Entry object.
    nmrStarExport = NmrStarExport(entry, nmrStarVersion = '3.1')
    nmrStarExport.createFile(outNmrStarFile)
    nmrStarExport.writeFile()
