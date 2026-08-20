# Exercise 1.3: create objects from top objects.
#
# - Create a new project object and top objects for the Affiliation and
#     Citation packages.
#
# - From the Affiliation top object, create a new person with a first name
#     and surname.
#
# - Also create a simple journal reference with a test title and journal
#     name from the Citation top object.
#
# - Then link the person to the citation - use the addAuthor() method
#     of the JournalCitation class.
#
# - Finally print the person's name and also the citation title starting
#     from the person object by finding the first citation for this person
#     with the relevant journal name.
#
# - Hint: use findFirstAuthorCitation to follow the person to citation link
#     in this case.
#

from memops.api.Implementation import MemopsRoot

if __name__ == '__main__':

  project = MemopsRoot(name = 'objTest')
  projectName = project.name

  # Create a new AffiliationStore top object.

  affStore = project.newAffiliationStore(name = projectName)

  # Create a new CitationStore top object.

  citStore = project.newCitationStore(name = projectName)

  # Create a new Person object from the AffiliationStore top object.

  person = affStore.newPerson(familyName = 'Smith', givenName = 'John')

  # Create a new JournalCitation object from the CitationStore top object.

  journal = citStore.newJournalCitation(title = 'My paper',
                                        journalAbbreviation = 'JMB')

  # Make a link between the person and the citation.

  journal.addAuthor(person)

  # Navigate from the person to his citation.

  person_citation = person.findFirstAuthorCitation(journalAbbreviation = 'JMB')

  # Print some info about the person and his citation.

  print 'Person: [%s]; Journal [%s]' % (person.familyName, person_citation.title)
