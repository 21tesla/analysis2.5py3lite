import os
import Tkinter

from memops.api.Implementation import (MemopsRoot, Url)
from memops.universal import Io as uniIo

from ccpnmr.format.converters.CnsFormat import CnsFormat

from ccpnmr.format.general.Util import (createExperiment,
                                        getRefExpFromOldExpType)

from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport

# Sub routine to show how to save a project to a user-defined location.
def saveProject(project, name = 'test', subDir = 'data'):

    # Find the backup repository and set the right path for it.
    backupRepos = project.findFirstRepository(name = 'backup')
    backupPath  = uniIo.joinPath(os.getcwd(), subDir, name + '_backup')

    # Create a URL for the backup repository.
    backupRepos.url = Url(path = backupPath)

    # Do the same for the repository containing the data in the project.
    projectRepos = project.findFirstRepository(name = 'userData')
    newPath      = uniIo.joinPath(os.getcwd(), subDir, name)

    # And also for the URL of this data.
    projectRepos.url = Url(path = newPath)

    # This saves the project to the path specified by the URL for the
    # 'userData' repository.
    project.saveModified()

# Sub routine to print out details of the repositories.
def printReposUrls(project):

    print '\nUrls:'

    for repos in project.sortedRepositories():
        print repos.name, repos.url.path

    print

if __name__ == '__main__':

    # This code follows on from ex7.py.
    project = MemopsRoot(name = 'constraintsTest')
    guiRoot = Tkinter.Tk()

    cnsObj = CnsFormat(project, guiRoot)

    cnsDir = '../data/cns'

    cnsCoord     = 'cns_1.pdb'
    cnsCoordFile = os.path.join(cnsDir, cnsCoord)

    # To make a molSystem object.
    chains = cnsObj.readSequence(cnsCoordFile)

    distConst     = 'n15noesy.tbl'
    distConstFile = os.path.join(cnsDir, distConst)

    # To make a structureGeneration object.
    constList = cnsObj.readDistanceConstraints(distConstFile)

    nmrConstStore = constList.nmrConstraintStore
    strucGen = nmrConstStore.findFirstStructureGeneration()

    cnsObj.linkResonances(forceDefaultChainMapping = 1,
                          globalStereoAssign = 1,
                          setSingleProchiral = 1,
                          setSinglePossEquiv = 1,
                          strucGen = strucGen)

    # As usual, create the BMRB Entry object.
    entryStore = project.newNmrEntryStore(name = project.name)
    entry = entryStore.newEntry(name = project.name)

    # Get the molSystem object, and link to it.
    molSystem = project.sortedMolSystems()[0]
    entry.setMolSystem(molSystem)

    # Link the structureGeneration object so that the restraints
    # are part of the BMRB Entry object.
    entry.addStructureGeneration(strucGen)

    # Make a test NOESY experiment.
    refExpType = 'noesy_hsqc_HNH.hhn'
    refExp = getRefExpFromOldExpType(project, refExpType)

    nmrExp = createExperiment(project, 'noesyTest', refExp)
    nmrExp.addMolSystem(molSystem)

    # Link the experiment to the BMRB Entry.
    entry.addExperiment(nmrExp)

    # Save the project to XML files using a sub routine - see above.
    saveProject(project, name = 'testSave', subDir = '../data/cnsOut')

    # Sub routine to print where data for the CCPN project is saved.
    printReposUrls(project)

    curDir     = os.path.abspath('../data')
    nmrStarDir = os.path.join(curDir, 'nmrStar')

    if not os.path.exists(nmrStarDir):
        os.mkdir(nmrStarDir)

    outNmrStarFile = os.path.join(nmrStarDir, 'nmrStar4.str')

    nmrStarExport = NmrStarExport(entry, nmrStarVersion = '3.1')
    nmrStarExport.createFile(outNmrStarFile)
    nmrStarExport.writeFile()
