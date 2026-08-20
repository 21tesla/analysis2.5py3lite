import os
import Tkinter

from memops.api.Implementation import MemopsRoot

from ccpnmr.format.converters.CnsFormat import CnsFormat

from ccpnmr.format.general.Util import (createExperiment,
                                        getRefExpFromOldExpType)

from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport

if __name__ == '__main__':

    # Similar start to code as that found in ex6.py.
    project = MemopsRoot(name = 'sampleTest')
    guiRoot = Tkinter.Tk()

    cnsObj = CnsFormat(project, guiRoot)

    nmrProject = project.currentNmrProject = \
        project.newNmrProject(name = project.name)
    strucGen = nmrProject.newStructureGeneration()

    cnsDir  = '../data/cns'
    cnsFile = os.path.join(cnsDir, 'cns_1.pdb')

    cnsObj.readCoordinates([cnsFile],
                           strucGen = strucGen,
                           minimalPrompts = 1,
                           linkAtoms = 0)

    # Find molecule made by readCoordinates and set isFinalised to True.
    molecule = project.findFirstMolecule()
    molecule.isFinalised = True

    # Create a classification store.
    classif = project.CurrentClassification = \
        project.newClassification(namingSystem = 'local')

    # Make a new sampleCategory object - required for sample objects.
    sampCat = classif.newSampleCategory(name = 'myCategory')

    # Also need a sampleStore object.
    sampStore = project.currentSampleStore = \
        project.newSampleStore(name = project.name)

    # And make the new sample from the sampleStore with the sampleCategory.
    keywds = {'name':             'mySample',
              'sampleCategories': [sampCat],
              'ph':               5.0}
    sample = sampStore.newSample(**keywds)

    # Reference store for making molecule components.
    refStore = project.currentRefSampleComponentStore = \
        project.newRefSampleComponentStore(name = project.name)

    # New molComponent object connected to the protein molecule - subclass
    # of AbstractComponent. Needed for making a new sample component.
    keywds2 = {'name':     molecule.name,
               'molecule': molecule,
               'molType':  'protein'}
    molComp = refStore.newMolComponent(**keywds2)

    # Concentration data for a sampleComponent for the protein molecule.
    keywds3 = {'concentration':      0.001,
               'concentrationError': 0.0002,
               'concentrationUnit':  'M'}
    sampComp = sample.newSampleComponent(refComponent = molComp, **keywds3)

    # Old way of labelling information for the protein molecule - see ex11b.py
    #   for the new way to do this (which will be seen in the NmrStar file).

    keywds4 = {'elementName':   'C',
               'incorporation': 0.99,
               'isUniform':     True,
               'labelName':     '13C',
               'labelType':     'isotope',
               'massNumber':    13}
    label13C = molComp.newLabel(**keywds4)

    keywds5 = {'elementName':   'N',
               'incorporation': 0.98,
               'isUniform':     True,
               'labelName':     '15N',
               'labelType':     'isotope',
               'massNumber':    15}
    label15N = molComp.newLabel(**keywds5)

    # New substance object for NaCl - subclass of AbstractComponent.
    # Note, a substance doesn't need a molecule unlike a molComponent.
    saltComp = refStore.newMolComponent(name    = 'NaCl',
                                        details = 'NaCl')

    # New sampleComponent object connected to the substance object.
    keywds6 = {'concentration':     0.1,
               'concentrationUnit': 'M'}
    sampComp2 = sample.newSampleComponent(refComponent = saltComp, **keywds6)

    # sampleConditionSet object - parent of individual sampleCondition objects.
    sampCondSet = nmrProject.newSampleConditionSet(name = 'sample_conditions')

    # Create two sampleCondition objects.
    keywds7 = {'value':     5.0,
               'condition': 'pH'}
    sampCondPh = sampCondSet.newSampleCondition(**keywds7)

    keywds8 = {'value':     298,
               'condition': 'Temperature',
               'unit':      'K'}
    sampCondTemp = sampCondSet.newSampleCondition(**keywds8)

    # Easy way to make a test NMR experiment for connecting to
    # sample information.
    refExpType = 'noesy_hsqc_HCNH.hhcn'
    refExp = getRefExpFromOldExpType(project, refExpType)
    nmrExp = createExperiment(project, '13c_15n_noesy_test', refExp)

    # Link the sample information to this experiment.
    nmrExp.setSample(sample)
    nmrExp.setSampleConditionSet(sampCondSet)

    # Find the molSystem to connect to the experiment.
    molSystem = project.findFirstMolSystem()
    nmrExp.addMolSystem(molSystem)

    # affiliationStore for NMR spectrometer manufacturer.
    affStore = project.currentAffilitationStore = \
        project.newAffiliationStore(name = project.name)

    manufacturer = affStore.newOrganisation(name = 'Varian')

    # instrumentStore for the NMR machine.
    instrStore = project.currentInstrumentStore = \
        project.newInstrumentStore(name = project.name)

    keywds9 = {'name':         'Varian_Inova_600',
               'manufacturer': manufacturer,
               'nominalFreq':  '600',
               'model':        'Inova',
               'protonFreq':   600.1}
    spec = instrStore.newNmrSpectrometer(**keywds9)

    # Connect the NMR experiment to the spectrometer.
    nmrExp.setSpectrometer(spec)

    # Add some methods and software - made from a parent methodStore object.
    methStore = project.currentMethodStore = \
        project.newMethodStore(name = project.name)

    keywds10 = {'name': 'assignment',
                'task': 'NMR peak assignment'}
    method1 = methStore.newMethod(**keywds10)

    keywds11 = {'name': 'peak picking',
                'task': 'NMR peak picking'}
    method2 = methStore.newMethod(**keywds11)

    keywds12 = {'name':       'Analysis',
                'vendorName': 'ccpNmr',
                'version':    '1.0',
                'methods':    [method1, method2]}
    software = methStore.newSoftware(**keywds12)

    # Make a new BMRB Entry.
    entryStore = project.currentNmrEntryStore = \
        project.newNmrEntryStore(name = project.name)
    entry = entryStore.newEntry(name = project.name)

    # Set the relevant links to the new BMRB Entry.
    entry.setMolSystem(molSystem)
    entry.addExperiment(nmrExp)

    # The plural of software is software - hence setSoftware requires
    # a set of software.
    entry.setSoftware([software])  # or entry.addSoftware(software)
    entry.addStructureGeneration(strucGen)

    curDir     = os.path.abspath('../data')
    nmrStarDir = os.path.join(curDir, 'nmrStar')

    if not os.path.exists(nmrStarDir):
        os.mkdir(nmrStarDir)

    outNmrStarFile = os.path.join(nmrStarDir, 'nmrStar3.str')

    nmrStarExport = NmrStarExport(entry, nmrStarVersion = '3.1')
    nmrStarExport.createFile(outNmrStarFile)
    nmrStarExport.writeFile()
