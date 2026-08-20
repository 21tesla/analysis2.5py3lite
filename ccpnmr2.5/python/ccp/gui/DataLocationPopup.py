
"""
======================COPYRIGHT/LICENSE START==========================

DataLocationPopup.py: <write function here>

Copyright (C) 2008 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
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

===========================REFERENCE END===============================
"""

from memops.gui.Util import createDismissHelpButtonList

from memops.editor.BasePopup import BasePopup

from ccp.gui.DataLocationFrame import DataLocationFrame

class DataLocationPopup(BasePopup):

  def __init__(self, parent, project, title = 'Edit Project Spectrum & Data Paths',
               help_msg = '', help_url = '', *args, **kw):

    self.help_msg = help_msg
    self.help_url = help_url
    self.waiting  = False
    BasePopup.__init__(self, parent, project, title=title, *args, **kw)

  def body(self, guiFrame):

    self.geometry('700x500')
    guiFrame.expandGrid(0,0)

    self.dataLocationFrame = DataLocationFrame(guiFrame, self.project,
                                               utilityButtons=True, grid=(0,0))

    buttonList = self.dataLocationFrame.buttonList
    buttonList.helpUrl = self.help_url
    buttonList.helpMsg = self.help_msg
    
    self.dataLocationFrame.administerNotifiers(self.registerNotify)

  def destroy(self):

    self.dataLocationFrame.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 2:
    print 'must specify project directory'
    sys.exit()

  path = sys.argv[1]

  import Tkinter
  from memops.gui.DataEntry import askDir, askFile
  from memops.gui.MessageReporter import showWarning
  from memops.general.Io import loadProject
  from memops.universal.Io import normalisePath

  # This needs to be above loadProject because that can pop up
  # dialogs, and it will create a root if one has not been created
  # already, and that will lead to the later code crashing
  r = Tkinter.Tk()

  path = normalisePath(path)
  askdir = lambda title, prompt, initial_value: askDir(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
  askfile = lambda title, prompt, initial_value: askFile(title, prompt,
              initial_value, parent=top, extra_dismiss_text='Skip')
  project = loadProject(path, showWarning=showWarning,
                        askDir=askdir, askFile=askfile)

  popup = DataLocationPopup(r, project)

  r.withdraw()
  r.mainloop()
