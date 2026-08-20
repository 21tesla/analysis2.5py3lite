
"""
======================COPYRIGHT/LICENSE START==========================

SaveProjectPopup.py: <write function here>

Copyright (C) 2008 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
import os
import Tkinter

from memops.universal.Io import joinPath, splitPath, normalisePath

from memops.general import Implementation, Io

"""
from memops.gui.Button import Button
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showError, showInfo, showYesNo
from memops.gui.Util import createDismissHelpButtonList
"""

from memops.editor.BasePopup import BasePopup
from memops.editor.SaveProjectFrame import SaveProjectFrame

defaultProjectFile = 'project.xml'

class SaveProjectPopup(BasePopup):

  def __init__(self, parent, project, title = 'Project : Save As', callback = None,
               help_msg = '', help_url = '', dismiss_text = '', *args, **kw):

    self.callback = callback
    self.help_msg = help_msg
    self.help_url = help_url
    self.dismiss_text = dismiss_text
    BasePopup.__init__(self, parent=parent, project=project, title=title, *args, **kw)

  def body(self, guiParent):

    self.geometry('600x600')

    frame = SaveProjectFrame(guiParent, self.project, self.callback,
              self.help_msg, self.help_url, self.dismiss_text, modal=self.modal)
    frame.grid(row=0, column=0, sticky=Tkinter.NSEW)
    self.save_frame = frame

  def __getattr__(self, name):

    try:
      return getattr(self.__dict__['save_frame'], name)
    except:
      raise AttributeError, "%s instance has no attribute '%s'" % (self.__class__.__name__, name)

