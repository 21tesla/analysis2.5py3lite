# Exercise 1.2: navigate a CCPN project.
#
# - Load up the previously made project using loadProject(<projectDir>)
#     in 'memops.general.Io'.
#
# - Look at the database top object you created using findFirst<Element>
#     with the parameter of the database name.
#
# - Get a sorted list of the taxonomy top objects and also find all NmrProjects.
#
# - Run the script in interactive mode using python -i <scriptname>.
#
# - Have a look at the database object in more detail using the dir() function.
#
# - Try to look at the output of some of the links and methods.
#
# - For the database object, follow the standard links (all objects have these)
#     to look at the:
#   - root
#   - memopsRoot
#   - parent
#   - topObject
#

from memops.general.Io import loadProject

if __name__ == '__main__':

  project = loadProject('topObjectTest')

  print 'Project: [%s]' % project, '\n'

  # Find the first database with a parameter.

  database = project.findFirstDatabase(name = 'UniProt')

  print 'Database: [%s]' % database, '\n'

  # Sorted taxonomy top objects.

  taxonomies = project.sortedTaxonomies()

  print 'Sorted taxonomies: [%s]' % taxonomies, '\n'

  # Find all NmrProjects.

  nmrProjects = project.findAllNmrProjects()

  print 'All Nmr Projects: [%s]' % nmrProjects, '\n'

  # The dir() function.

  print 'Database dir:', dir(database), '\n'

  # root
  print 'Database root:', database.root, '\n'

  # memopsRoot
  print 'Database memopsRoot:', database.memopsRoot, '\n'

  # parent
  print 'Database parent:', database.parent, '\n'

  # topObject
  print 'Database top object:', database.topObject, '\n'
