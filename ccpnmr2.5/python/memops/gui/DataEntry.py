
"""
======================COPYRIGHT/LICENSE START==========================

DataEntry.py: <write function here>

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
import memops.gui.QueryDialogBox as QueryDialogBox


from memops.gui.FileSelectPopup import FileSelectPopup

def askPassword(title, prompt, parent = None):

  return QueryDialogBox.askPassword(title, prompt, parent=parent)

def askString(title, prompt, initial_value = '', parent = None):

  return QueryDialogBox.askString(title, prompt,initialvalue=initial_value,
                                  parent=parent)


def askInteger(title, prompt, initial_value = '', min_value = None,
               max_value = None, parent = None):

  return QueryDialogBox.askInteger(title, prompt, initialvalue=initial_value,
                                   minvalue=min_value, maxvalue=max_value, parent=parent)

def askFloat(title, prompt, initial_value = '', min_value = None,
             max_value = None, parent = None):

  return QueryDialogBox.askFloat(title, prompt, initialvalue=initial_value,
                                  minvalue=min_value, maxvalue=max_value, parent=parent)


def askFile(title, prompt, initial_value = '', parent = None,
            dismiss_text='Cancel', extra_dismiss_text = ''):

  if (parent):
    popup = FileSelectPopup(parent, title=title, prompt=prompt, show_file=True,
                            dismiss_text=dismiss_text,
                            extra_dismiss_text=extra_dismiss_text,
                            file=initial_value)
      
    file = popup.getFile()
    popup.destroy()
    return file
  else:
    return askString(title, prompt, initial_value)

def askDir(title, prompt, initial_value = '', parent = None,
           dismiss_text='Cancel', extra_dismiss_text = '', default_dir = None):

  if (parent):
    popup = FileSelectPopup(parent, title=title, prompt=prompt, show_file=False,
                            dismiss_text=dismiss_text,
                            extra_dismiss_text=extra_dismiss_text,
                            file=initial_value, default_dir = default_dir)

    dir = popup.getDirectory()
    popup.destroy()
    return dir
  else:
    return askString(title, prompt, initial_value)

class DataEntry:

  def askPassword(self, title, prompt, initial_value = '', parent = None, *args, **kw):
 
    return askPassword(title, prompt, initial_value, parent)
 
  def askString(self, title, prompt, initial_value = '', parent = None, *args, **kw):
 
    return askString(title, prompt, initial_value, parent)
 
  def askInteger(self, title, prompt, initial_value = '', min_value = None,
                 max_value = None, parent = None, *args, **kw):
 
    return askInteger(title, prompt, initial_value, min_value, max_value, parent)
 
  def askFloat(self, title, prompt, initial_value = '', min_value = None,
               max_value = None, parent = None, *args, **kw):
 
    return askFloat(title, prompt, initial_value, min_value, max_value, parent)
 
  def askFile(self, title, prompt, initial_value = '', parent = None,
              dismiss_text='Cancel', extra_dismiss_text = '', *args, **kw):
 
    return askFile(title, prompt, initial_value, parent)
 
  def askDir(self, title, prompt, initial_value = '', parent = None,
             dismiss_text='Cancel', extra_dismiss_text = '', default_dir = None, *args, **kw):
 
    return askDir(title, prompt, initial_value, parent, default_dir = default_dir)
 
dataEntry = DataEntry()

if (__name__ == '__main__'):
 
  import Tkinter

  r = Tkinter.Tk()
  print dataEntry.askString('ask string title', 'ask string prompt')
  print dataEntry.askInteger('ask integer title', 'ask integer prompt')
  print dataEntry.askFloat('ask float title', 'ask float prompt')
  print dataEntry.askFile('ask file title', 'ask file prompt', parent=r)
  print dataEntry.askDir('ask dir title', 'ask dir prompt', parent=r)
