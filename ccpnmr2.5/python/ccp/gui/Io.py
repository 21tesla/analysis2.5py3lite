"""
======================COPYRIGHT/LICENSE START==========================

Io.py: Allows loading of projects.

Copyright (C) 2008 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""

import os

from memops.universal import Io as uniIo

from memops.general import Io as genIo

from memops.gui.DataEntry import askDir, askFile
from memops.gui.MessageReporter import showWarning

from ccp.gui.DataLocationPopup import DataLocationPopup

from ccp.general import Io as ccpGenIo

from ccp.api.general.DataLocation import NumericMatrix
from ccp.api.general.DataLocation import MimeTypeDataStore

def loadProject(parent, path, projectName=None):

  path = uniIo.normalisePath(path)
  askdir = lambda title, prompt, initial_value: askDir(title, prompt,
              initial_value, parent=parent, extra_dismiss_text='Skip')
  askfile = lambda title, prompt, initial_value: askFile(title, prompt,
              initial_value, parent=parent, extra_dismiss_text='Skip')
  project = genIo.loadProject(path, showWarning=showWarning,
                              askDir=askdir, askFile=askfile)

  # now check dataStores
  # delete those that are not used
  # and otherwise check path to see if exists

  dataStores = []
  for dataLocationStore in project.dataLocationStores:
    for dataStore in dataLocationStore.dataStores:
      if isinstance(dataStore, NumericMatrix) and not dataStore.nmrDataSources:
        print 'deleting dataStore %s with path %s' % (dataStore, dataStore.fullPath)
        dataStore.delete()
      elif isinstance(dataStore, MimeTypeDataStore) and not dataStore.nmrDataSourceImages:
        print 'deleting dataStore %s with path %s' % (dataStore, dataStore.fullPath)
        dataStore.delete()
      else:
        dataStores.append(dataStore)

  badDataStores = [ dataStore for dataStore in dataStores if not os.path.exists(dataStore.fullPath) ]

  if badDataStores:
    # find DataUrls involved 
    dataUrls = set(dataStore.dataUrl for dataStore in badDataStores)
    startDir = project.packageLocator.findFirstRepository().url.dataLocation
    
    for dataUrl in dataUrls:
      if not dataUrl.dataStores.difference(badDataStores):
        # all DataStores for this DataUrl are bad
        # we can make changes without affecting 'good' DataStores
 
        # Look for an obvious place the data may have moved to
        dataStores =  dataUrl.sortedDataStores()
        fullPaths = [dataStore.fullPath for dataStore in dataStores]
        baseDir, newPaths = uniIo.suggestFileLocations(fullPaths,
                                                       startDir=startDir)
 
        if baseDir is not None:
          # We have a file location that fits all missing files.
          # Change dataStores to use it
          print 'WARNING, resetting data locations to: \n%s\n' % baseDir
 
          ccpGenIo.changeDataStoreUrl(dataStores[0], baseDir)
          for ii,dataStore in enumerate(dataStores):
            dataStore.path = newPaths[ii]
    
    if [dataStore for dataStore in dataStores if not os.path.exists(dataStore.fullPath)]:
      popup = DataLocationPopup(parent, project, modal=True)
      popup.destroy()

  return project
