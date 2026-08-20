
"""
======================COPYRIGHT/LICENSE START==========================

Util.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
import sys, os
import types

import Tkinter

from memops.gui.Base import getPopup
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.WebBrowser import WebBrowser
import memops.gui.HelpPopup

from memops.universal.Io import getTopDirectory

gfxDir = os.path.join(getTopDirectory(),'python','memops','gui','graphics')

# creates quit button

def createQuitButton(parent, text=''):

  if not text:
    text = ''

  b = Button(parent, text=text, command=sys.exit)
  
  if not text:
    b.cancelIcon = Tkinter.PhotoImage(file=os.path.join(gfxDir,'cancel.gif'))
    b.config(image=b.cancelIcon, activebackground=b.cget('bg'))

  return b

# creates dismiss button
# if dismiss_cmd is not set then assumes that parent is in a BasePopup

def createDismissButton(parent, dismiss_text = '', dismiss_cmd = None):

  if (not dismiss_text):
    dismiss_text = ''

  if (not dismiss_cmd):
    popup = parent
    while (not hasattr(popup, 'top')):
      popup = popup.parent
    dismiss_cmd = popup.close

  b = Button(parent, text=dismiss_text, command=dismiss_cmd)
  if not dismiss_text:
    b.cancelIcon = Tkinter.PhotoImage(file=os.path.join(gfxDir,'cancel.gif'))
    b.config(image=b.cancelIcon, activebackground=b.cget('bg'))

  return b

# creates dismiss and help buttons and other ones (in front) if desired
# if dismiss_cmd is not set then assumes that parent is in a BasePopup

def createDismissHelpButtonList(parent, texts = None, commands = None,
                          direction=Tkinter.HORIZONTAL,
                          dismiss_text = '', help_text = '',
                          help_msg = '', help_url = '',
                          buttonBorderwidth=True, expands=True, 
                          dismiss_cmd = None,
                          webBrowser = None,
                          tipTexts = None,
                          *args, **kw):
  if texts is None:
    texts = []

  if commands is None:
    commands = []

  if not dismiss_text:
    dismiss_text = ''
 
  if not help_text:
    help_text = ''

  if not tipTexts:
    tipTexts = len(texts) * [None]

  texts = list(texts) + [dismiss_text, help_text]
  tipTexts = list(tipTexts) + [dismiss_text, help_text]

  popup = getPopup(parent)

  if not webBrowser:
    webBrowser = WebBrowser(popup.top, popup=popup)

  if (not dismiss_cmd):
    dismiss_cmd = popup.close

  if (help_url):
    help_cmd   = lambda url=help_url: webBrowser.open(url)
  else:
    help_cmd = lambda top=popup.top, message=help_msg: memops.gui.HelpPopup.showHelpText(top, message, popup=popup)

  if type(commands) is types.DictType:
    commands = commands.copy()
    commands[dismiss_text] = dismiss_cmd
    commands[help_text] = help_cmd
  else:
    commands = list(commands) + [dismiss_cmd, help_cmd]

  button_list = ButtonList(parent, texts=texts, tipTexts=tipTexts,
                           commands=commands, buttonBorderwidth=buttonBorderwidth,
                           expands=expands, direction=direction, *args, **kw)
  
  
  if not dismiss_text:
    button_list.cancelIcon = Tkinter.PhotoImage(file=os.path.join(gfxDir,'cancel.gif'))
    button_list.buttons[-2].config(image=button_list.cancelIcon, activebackground=button_list.cget('bg'))
  
  if not help_text:
    button_list.helpIcon   = Tkinter.PhotoImage(file=os.path.join(gfxDir,'help.gif'))
    button_list.buttons[-1].config(image=button_list.helpIcon, activebackground=button_list.cget('bg'))
    

  return button_list

#
# Same as createDismissHelpButtonList, but with dismiss button taken out
#

def createHelpButtonList(parent, texts = None, commands = None,
                          direction=Tkinter.HORIZONTAL,
                          help_text = '',
                          help_msg = '', help_url = '',
                          buttonBorderwidth=True, expands=True,
                          webBrowser = None,
                          *args, **kw):
  
  if (texts is None):
    texts = []

  if (commands is None):
    commands = []
 
  if not help_text:
    help_text = ''

  texts = list(texts) + [help_text]

  popup = getPopup(parent)

  if not webBrowser:
    webBrowser = WebBrowser(popup.top, popup=popup)

  if (help_url):
    help_cmd   = lambda url=help_url: webBrowser.open(url)
  else:
    help_cmd = lambda top=popup.top, message=help_msg: memops.gui.HelpPopup.showHelpText(top, message, popup=popup)

  if (type(commands) is types.DictType):
    commands = commands.copy()
    commands[help_text] = help_cmd
  else:
    commands = list(commands) + [help_cmd]

  button_list = ButtonList(parent, texts=texts, commands=commands, buttonBorderwidth=buttonBorderwidth,
                           expands=expands, direction=direction, *args, **kw)

  
  if not help_text:
    button_list.helpIcon = Tkinter.PhotoImage(file=os.path.join(gfxDir,'help.gif'))
    button_list.buttons[-1].config(image=button_list.helpIcon)

  return button_list

