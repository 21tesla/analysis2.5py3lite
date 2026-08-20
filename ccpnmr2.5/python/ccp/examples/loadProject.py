"""
##################
# Example script #
##################

Shows how to load a CCPN project from an XML project file. Also resets the path names.

Level: basic

Contact: Wim Vranken <wim@ebi.ac.uk>
"""

#
# Standard Python import
#

import os

#
# Import the Implementation package for creating Urls
#

from memops.api import Implementation
from memops.universal import Io as uniIo

#
# Import the CCPN xml file loader
#

from memops.general.Io import loadProject

#
# Function to get the CCPN top directory
#

from memops.universal.Io import getTopDirectory

#
# Function to join directories
#
# NBNB TBD currently this must be sourced from the directory taht contains it


from memops.universal.Io import joinPath

###################
# Main of program #
###################

if __name__ == "__main__":

  #
  # Load a CCPN project
  #
  
  projectName = 'test'
  projectDir = os.path.join(os.getcwd(),'local',projectName)

  ccpnProject = loadProject(projectDir,projectName)

  #
  # Because the original project was saved in a different directory structure,
  # the paths have to be reset first... here I am resetting everything 'hardcoded'
  # based on the examples and/or data directory.
  #

  localDataPaths = {}
  localDataPaths['examples'] = os.path.abspath('.')
  localDataPaths['data'] = joinPath(getTopDirectory(),'data')

  for repository in ccpnProject.repositories:
    for findName in ('examples','data'):
      findNameIndex = repository.url.path.find(findName)
      if  findNameIndex>= 0:
        repository.url = Implementation.Url(path = uniIo.joinPath(localDataPaths[findName],repository.url.path[findNameIndex + len(findName) + 1:]))
        break

  #
  # Print out some objects from this project
  #
  
  print "Project object:", ccpnProject
  print "Project name:", ccpnProject.name
  
  #
  # Note that there is a pause after the above is printed - this is because of
  # 'lazy loading' the API will now have to read the MolSystem.xml file before it
  # can print the information in there.
  #
  for mol in ccpnProject.molecules:
    print mol.name
    print mol.sortedMolResidues()
    
  print "List of molSystem objects:", ccpnProject.sortedMolSystems()
  for molSys in ccpnProject.sortedMolSystems():
    print "MolSystem name:", molSys.name
    print "MolSystem code:", molSys.code
  
  print "List of molecule objects:", ccpnProject.sortedMolecules()
  for molecule in ccpnProject.sortedMolecules():
    print "Molecule name:", molecule.name
