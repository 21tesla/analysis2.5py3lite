# Exercise 1.4: load an old CCPN project.
#
# - Load up the project from exercise 1.1 .
#
# - Find the 'userData' repository for this project and save the
#     project to a new location on the disk.
#
# - Hint: specify the new location as the Url link object for this repository.
#
# - Optional: create a new repository with Url path '../data'.
#
# - Then add this repository to the PackageLocator object with
#     targetName = 'ccp.molecule.ChemComp'.
#
# - There are some carbohydrate chemcomps in '../data/ccp/molecule/ChemComp'.
#

# Import the standard python module called os.

import os

# Import the loadProject() function.

from memops.general.Io import loadProject
from memops.universal import Io as uniIo

# Import the Url class.

from memops.api.Implementation import Url

if __name__ == '__main__':

  # Load the old project - need to specify the correct directory name.

  project = loadProject('topObjectTest')

  print 'Project: [%s]\n' % project

  # Find the right repository that contains the project data.

  projectRepos = project.findFirstRepository(name = 'userData')

  # Specify a new location to save the project in.

  newPath = uniIo.normalisePath(os.path.join(os.path.abspath('.'), 
                                             'newDir', 'newSubDir'))

  # Make a new Url pointing to the new location and link it to the
  # 'userData' repository.

  projectRepos.url = Url(path = newPath)

  # Print the location of all project repositories.

  print 'Repositories:\n%s' % (
    '\n'.join(repos.name + ': ' + repos.url.path
              for repos in project.sortedRepositories() ) )

  carbDataPath = uniIo.normalisePath(os.path.join(os.path.abspath('..'), 'data'))

  print '\nLocation of local carbohydrate ChemComps: [%s]' % carbDataPath

  carbDataUrl = Url(path = carbDataPath)

  carbDataRepos = project.newRepository(name = 'carbData',
                                        url = carbDataUrl)

  chemPackLoc = project.findFirstPackageLocator(targetName = 'ccp.molecule.ChemComp')

  chemPackLoc.addRepository(carbDataRepos)

  # Save the project to the new location.

  #project.saveAll()
