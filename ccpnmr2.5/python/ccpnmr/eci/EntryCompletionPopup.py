
"""
======================COPYRIGHT/LICENSE START==========================

EntryCompletionPopup.py: Part of the CcpNmr Analysis program

Copyright (C) 2005 Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: tjs23@cam.ac.uk
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

from memops.gui.ButtonList import UtilityButtonList
from memops.editor.BasePopup import BasePopup

from ccpnmr.eci.EntryCompletionFrame import EntryCompletionFrame

def testEditNmrEntryPopup(argServer):

  popup = EntryCompletionPopup(argServer.parent, argServer.getProject())
  popup.open()


class EntryCompletionPopup(BasePopup):

  """
  **Setup Deposition of Coordinate and NMR Data to PDB and BMRB**
  
  The popup version of CcpNmr Entry Completion Interface (ECI), available from
  Analysis in the Structure menu.  With it, you can complete all the necessary
  information required for PDB and BMRB data depositions by adding an "Entry"
  object to your CCPN project. An "Entry" object contains all the information that
  you wish to deposit with your submission. You can also select chemical shift
  lists, peak lists, structural restraints, ensembles, etc. and all at the click
  of a button. In addition, you can add all the meta data that is required for
  submissions to the PDB and BMRB. This can be done securely on your desktop
  computer over the duration of your NMR project.

  For more documentation see the `ECI section of the PDBe pages at the EBI web
  site <http://www.ebi.ac.uk/pdbe/docs/pdbe_nmr_deposition/eci.html>`_.
  """

  def __init__(self, parent, project, *args, **kw):
  
    self.project = project             
    BasePopup.__init__(self, parent, title='Deposition Entry Completion Interface', **kw)
    
    
  def body(self, guiFrame):

    self.geometry('950x700')
    
    guiFrame.expandGrid(0,0)

    frame = EntryCompletionFrame(guiFrame, basePopup=self, grid=(0,0))
    frame.updateAll()

    utilButtons = UtilityButtonList(frame.tabbedFrame.sideFrame)
    utilButtons.grid(row=0, column=0, sticky='e')

