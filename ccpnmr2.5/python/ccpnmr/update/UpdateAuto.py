"""
======================COPYRIGHT/LICENSE START==========================

UpdateAuto.py: Part of the CcpNmr Update program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

===========================REFERENCE END===============================

"""
from ccpnmr.update.UpdateAgent import UpdateAgent, UPDATE_SERVER_LOCATION, UPDATE_DIRECTORY, UPDATE_DATABASE_FILE

def updateAuto(serverLocation=UPDATE_SERVER_LOCATION,
               serverDirectory=UPDATE_DIRECTORY,
               dataFile=UPDATE_DATABASE_FILE):

  updateAgent = UpdateAgent(serverLocation, serverDirectory, dataFile, isGraphical=False)
  server = updateAgent.server
  if server:
      server.getFileUpdates()
      if not server.fileUpdates:
        server.getFileUpdates()
      if server.fileUpdates:
        updateAgent.installNewUpdates()  
      else:
        print 'No updates'
    
if __name__ == '__main__':
  import os

  updateAuto()

  # NOTE:ED - seems to need this to clean up in Windows10
  os._exit(0)