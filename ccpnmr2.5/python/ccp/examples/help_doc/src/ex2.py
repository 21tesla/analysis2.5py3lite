from memops.api.Implementation import MemopsRoot

# No longer need the following package imports as can use methods
# in parent objects to create child objects like NmrEntryStore or Entry.
#from ccp.api.nmr.NmrEntry import (NmrEntryStore, Entry)

if __name__ == '__main__':

    # Still need to create the CCPN project from an import.
    project = MemopsRoot(name = 'entryTest2')

    # Can use methods to get object attributes.
    projectName = project.getName()
    print 'Project [%s]' % projectName

    # New way to create an NmrEntryStore instance.
    entryStore = project.newNmrEntryStore(name = projectName)
    entryStoreName = entryStore.getName()
    print 'Entry store [%s]' % entryStoreName

    # New way to create an Entry instance.
    keywds = {'title':          'test title',
              'bmrbProcessing': 'test submission',
              'details':        'test details'}
    entry = entryStore.newEntry(name = entryStoreName + 'A', **keywds)

    entryName     = entry.getName()
    entryTitle    = entry.getTitle()
    entryBmrbProc = entry.getBmrbProcessing()
    entryDetails  = entry.getDetails()

    print 'Entry [%s] [%s] [%s] [%s]\n' % (entryName, entryTitle,
                                           entryBmrbProc, entryDetails)

    # Make a second Entry with slightly different attributes.
    keywds2 = {'title':          'test title2',
               'bmrbProcessing': 'test submission',
               'details':        'test details'}
    entry2 = entryStore.newEntry(name = entryStoreName + 'B', **keywds2)

    # Entries are stored as a frozenset in NmrEntryStore.
    # Can use various methods to retrieve entries from this store.
    firstEntry = entryStore.findFirstEntry(name = entryStoreName + 'B')
    print 'First Entry with name [%s] is [%s]\n' % (entryStoreName + 'B',
                                                    firstEntry)

    findAllEntries = entryStore.findAllEntries(details = entryDetails)
    print "All entries with details [%s] as a 'set':\n[%s]\n" % (entryDetails,
                                                                 findAllEntries)

    allEntries = entryStore.getEntries()
    print "'frozenset' of all entries in a random order:\n[%s]\n" % allEntries

    sortedEntries = entryStore.sortedEntries()
    print "'list' of all entries sorted:\n[%s]\n" % sortedEntries
