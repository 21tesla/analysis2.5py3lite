# Package imports from CCPN API.
# These must be in the global scope.
from memops.api.Implementation import MemopsRoot
from ccp.api.nmr.NmrEntry import (NmrEntryStore, Entry)

if __name__ == '__main__':

    # Create a CCPN project (called a MemopsRoot class in the CCPN API).
    # Must specify the attribute 'name' when created.
    project = MemopsRoot(name = 'entryTest')

    # Non-API way to get object attributes.
    projectName = project.name
    print 'Project [%s]' % projectName

    # Create an NmrEntryStore that can hold multiple BMRB Entry objects.
    entryStore = NmrEntryStore(project, 
                               name = projectName)

    entryStoreName = entryStore.name
    print 'Entry store [%s]' % entryStoreName
    print 'Current store: [%s]' % project.currentNmrEntryStore

    # Create a BMRB Entry object, which can hold all the data for a BMRB
    # deposition.
    entry = Entry(entryStore, name = entryStoreName)
    entryName = entry.name
    print 'Entry [%s]\n' % entryName
