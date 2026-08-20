# Exercise 1.1: make a CCPN project and some associated top objects.
#
# - Create a new MemopsRoot (project) object.
#
# - From this, create top objects for the following packages:
#   - NmrEntry
#   - Nmr
#   - Affiliation
#   - Citation
#   - DbRef
#   - Taxonomy
#   - Classification
#
# - Use the new<Class> factory methods available to the MemopsRoot class.
#
# - Print out the objects and their names for those that you create.
#
# - Note the package locations for all the objects.
#
# - Save the project to XML files on disk and have a quick look at some
#     of the files that are created - use the method saveModified() .
#
# - Try and run the script twice and see what happens.
#
# - Hint: Check to see what the mandatory attributes are for these top objects.
#

# Import the MemopsRoot class from the API.

from memops.api.Implementation import MemopsRoot

if __name__ == '__main__':

  # New CCPN project - called a MemopsRoot.

  project = MemopsRoot(name = 'topObjectTest')

  # Use the new<Class> factory methods available to the MemopsRoot class.

  # NmrEntry
  entryStore = project.newNmrEntryStore(name = 'myEntry')

  # Nmr
  nmrProject = project.newNmrProject(name = 'myNmrProject')

  # Affiliation
  affStore = project.newAffiliationStore(name = 'myAffiliation')

  # Citation
  citStore = project.newCitationStore(name = 'myCitation')

  # DbRef
  database = project.newDatabase(name = 'UniProt')

  # Taxonomy
  taxonomy = project.newTaxonomy(name = 'myTaxonomy')

  # Classification
  classification = project.newClassification(namingSystem = 'myClassification')

  # Print the actual objects as strings.

  print '''
Project: [%s]; \nEntryStore: [%s]; \nNmrProject: [%s]\nAffiliation: [%s]; 
Citation: [%s]; \nDatabase: [%s]; \nTaxonomy: [%s]; \nClassification: [%s].
''' % (project, entryStore, nmrProject, affStore, citStore, database, taxonomy, classification)

  # Print the object names.

  print '''
Project: [%s]; \nEntryStore: [%s]; \nNmrProject: [%s]; \nAffiliation: [%s]; 
Citation: [%s]; \nDatabase: [%s]; \nTaxonomy: [%s]; \nClassification: [%s].
''' % (project.name, entryStore.name, nmrProject.name, affStore.name,
       citStore.name, database.name, taxonomy.name, classification.namingSystem)

  # Save the project.

  project.saveModified()
