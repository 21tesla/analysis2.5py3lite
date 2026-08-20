import os

# This contains the application data class.
from memops.api import Implementation
from memops.api.Implementation import MemopsRoot

from ccpnmr.format.converters.FastaFormat import FastaFormat

from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport

if __name__ == '__main__':

    project = MemopsRoot(name = 'alignTest')

    # Create a BMRB Entry for export to the NMR-STAR file.
    entryStore = project.newNmrEntryStore(name = project.name)
    nmrEntry   = entryStore.newEntry(name = project.name)

    # Application data objects allow you to store and associate extra data
    # in the CCPN project that does not fit into the rest of the data model.

    # Format of an application data object has the application that the data
    # is derived from (e.g., XEasy, NmrStar); a keyword to describe the type
    # of data that is present (e,g., sequence, chemical shift) and the value
    # of the data.

    # For this example, we are reading a Fasta file.
    seqDir  = '../data/seq'
    seqFile = os.path.join(seqDir, 'fasta.seq')
    
    keywds = {'application': 'Fasta format',
              'keyword':     'sequence file',
              'value':       seqFile}

    # Make and add this application data to the BMRB Entry.
    appData = Implementation.AppDataString(**keywds)
    nmrEntry.addApplicationData(appData)

    # Print it.
    entryAppData = nmrEntry.findFirstApplicationData()

    print 'Application data: [%s]' % entryAppData

    # Find it based on the application and keyword values.
    keywds2 = {'application': 'Fasta format',
               'keyword':     'sequence file'}

    seqFileFromAppData = nmrEntry.findFirstApplicationData(**keywds2).value

    # Read in a sequence to create a molSystem object - can also set the
    # Tkinter object to None for using all text-based format converter sub
    # routines.
    fastaObj = FastaFormat(project, None)

    fastaObj.readSequence(seqFileFromAppData, minimalPrompts = 1)

    # Connect the molSystem to the BMRB Entry.
    nmrEntry.setMolSystem(project.sortedMolSystems()[0])

    # Also create a new study object.
    keywordList = ['NMR', 'Structure determination']
    study = entryStore.newStudy(name = project.name,
                                keywords = keywordList)

    # And connect to the Entry object.
    nmrEntry.setStudy(study)

    affStore = project.newAffiliationStore(name = project.name)

    orgn  = affStore.newOrganisation(name = 'myDepartment')
    group = orgn.newGroup(name = 'myGroup')

    keywds = {'familyName': 'Smith',
              'givenName':  'John'}
    person = affStore.newPerson(**keywds)

    personInGroup = person.newPersonInGroup(group = group)

    nmrEntry.addContactPerson(person)
    nmrEntry.addAuthor(person)
    nmrEntry.addLaboratory(group)

    citStore = project.newCitationStore(name = project.name)

    keywds2 = {'title':               'myPaper',
               'journalAbbreviation': 'JMB'}
    citation = citStore.newJournalCitation(**keywds2)

    citation.addAuthor(person)
    nmrEntry.setPrimaryCitation(citation)

    # Create a new database entry - the top hit in UniProt from a
    # Blast search of the previously read in sequence.
    dbStore = project.newDatabase(name='UniProt')
    dbEntry = dbStore.newEntry(code = 'P27927')

    # Find the molecule object in the molSystem.
    molecule = nmrEntry.molSystem.sortedChains()[0].molecule

    # Some of the results of the Blast search for a sequence alignment object.
    keywds3 = {'alignLength':      360,
               'alignmentProgram': 'NCBI BlastP 2.2.17',
               'alignmentScore':   1507,
               'dbRefAlignBegin':  231,
               'dbRefAlignEnd':    586,
               'homologyRatio':    0.79,  # 79% sequence identity.
               'nIdentical':       287,
               'nPositive':        304}

    # Add a new alignment object to this molecule.
    molecule.newAlignment(dbRef = dbEntry, **keywds3)

    # Reference store for making molecule components.
    refSampCompStore = project.newRefSampleComponentStore(name = project.name)

    # New molecule component object for a molecule of type protein.
    keywds4 = {'name':    'testComp',
               'molType': 'protein'}
    molComp = refSampCompStore.newMolComponent(**keywds4)

    # More of the Blast results to link to the molecule component object.
    keywds5 = {'alignBegin':      1,
               'alignEnd':        338,
               'dbRefAlignBegin': 231,
               'dbRefAlignEnd':   586}

    # Link the molecule component to this UniProt cross reference.
    compDbRef = molComp.newComponentDbRef(dbRef = dbEntry, **keywds5)

    # Get the molResidues that begin and end the alignment.
    startResidue = molecule.findFirstMolResidue(seqCode = compDbRef.alignBegin)
    endResidue   = molecule.findFirstMolResidue(seqCode = compDbRef.alignEnd)

    # Create the molSeqFragment object for this alignment.
    limitResidues = [startResidue, endResidue]

    keywds6 = {'alignments':    molecule.sortedAlignments(),
               'limitResidues': limitResidues}
    molSeqFrag = molecule.newMolSeqFragment(**keywds6)

    # The molecule must be finalised to add molecule components.
    molecule.isFinalised = True

    # Add the previously made molComp object to the molecule.
    molecule.addMolComponent(molComp)

    # Taxonomy storage class.
    tax = project.newTaxonomy(name = project.name)

    # Create a natural source from this taxonomy object.
    keywds7 = {'organismName':   'Maize',
               'scientificName': 'Zea mays'}
    natSource = tax.newNaturalSource(**keywds7)

    # Set some links to the natural source.

    # Entry objects are linked to an entryMolecule object.
    #   - The entryMolecule is linked to a molecule
    #     that is also linked to a naturalSource object.
    #   - The entryMolecule is also linked directly to an
    #     experimentalSource object, and a productionMethod
    #     attribute.

    molecule.setNaturalSource(natSource)

    keywds8 = {'molecule':         molecule,
               'productionMethod': 'recombinant technology'}

    entryMol = nmrEntry.newEntryMolecule(**keywds8)

    molComp.setNaturalSource(natSource)
    molSeqFrag.setNaturalSource(natSource)

    # Also create an experimental source and link to the
    #   entryMolecule object.
    keywds9 = {'organismName':   'E. coli',
               'scientificName': 'Escherichia coli'}
    expSource = tax.newNaturalSource(**keywds9)

    entryMol.setExperimentalSource(expSource)

    curDir     = os.path.abspath('../data')
    nmrStarDir = os.path.join(curDir, 'nmrStar')

    if not os.path.exists(nmrStarDir):
        os.mkdir(nmrStarDir)

    outNmrStarFile = os.path.join(nmrStarDir, 'nmrStar2.str')

    nmrStarExport = NmrStarExport(nmrEntry, nmrStarVersion = '3.1')
    nmrStarExport.createFile(outNmrStarFile)
    nmrStarExport.writeFile()
