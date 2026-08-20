from memops.api.Implementation import MemopsRoot

if __name__ == '__main__':

    # CCPN project.
    project = MemopsRoot(name = 'personTest')

    # Create a BMRB Entry.
    entryStore = project.currentNmrEntryStore = \
        project.newNmrEntryStore(name = project.name)
    entry = entryStore.newEntry(name = project.name)

    # Create an affiliationStore and alias for it.
    affStore = project.currentAffiliationStore = \
        project.newAffiliationStore(name = project.name)

    # Create an organisation and group.
    orgn  = affStore.newOrganisation(name = 'myDepartment')
    group = orgn.newGroup(name = 'myGroup')

    # Create a new person.
    keywds = {'familyName': 'Smith',
              'givenName':  'John'}
    person = affStore.newPerson(**keywds)

    # Set the group for this person. PersonInGroup is an intermediate
    # class to avoid many-to-many relationships.
    personInGroup = person.newPersonInGroup(group = group)

    # Print some basic information.
    personName = person.givenName + ' ' + person.familyName
    orgName    = person.findFirstPersonInGroup().group.organisation.name
    print 'Person: [%s] is in Organisation [%s]' % (personName, orgName)

    # Set the person and group in the BMRB Entry.
    entry.addContactPerson(person)
    entry.addAuthor(person)
    entry.addLaboratory(group)

    # Create a citationStore and alias for it.
    citStore = project.currentCitationStore = \
        project.newCitationStore(name = project.name)

    # Create a citation.
    keywds2 = {'title': 'myPaper',
               'journalAbbreviation': 'JMB'}
    citation = citStore.newJournalCitation(**keywds2)

    # Use the abstract super class Citation to find the JournalCitation.
    myCit = citStore.findFirstCitation(className = 'JournalCitation')

    # Set some links for this citation to other objects in the project.
    myCit.addAuthor(person)

    # Also, see what happens if you comment out the following line.
    entry.setPrimaryCitation(myCit)

    # Print some more information.
    entryName = entry.name
    primCit   = entry.primaryCitation  # Shortcut to the previously made link.
    citAbbrev = primCit.journalAbbreviation

    # Here, the authors link is an ordered tuple and not a set.
    authorName = primCit.authors[0].familyName
    print 'Entry [%s] has a primary citation in [%s] authored by [%s]\n' % (
        entryName, citAbbrev, authorName)
