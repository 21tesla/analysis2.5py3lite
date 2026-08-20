
"""
======================COPYRIGHT/LICENSE START==========================

KeysymList.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

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

import string

from memops.gui.PulldownList import PulldownList

no_keysym = 'None'

class KeysymList(PulldownList):
 
  def __init__(self, parent, callback = None, *args, **kw):
 
    fkey_entries = [ 'F%d' % (n+1) for n in range(12) ]

    misc_entries = [ \
      'Home',
      'End',
      'Prior',
      'Next',
      'Up',
      'Down',
      'Right',
      'Left',
      'Delete',
    ]

    # remove 0-9 because used for spectrum shortcuts
    texts      = [no_keysym]
    categories = [None,]
    objects    = [None,]
    for category, sublabels in (('a-z', string.lowercase[:26]),
                               ('A-Z', string.uppercase[:26]),
                               #('0-9', string.digits),
                               ('F keys', fkey_entries),
                               ('Misc', misc_entries)):
      
      for text in sublabels:
        texts.append(text)
        objects.append(text)
        categories.append(category)
      

    PulldownList.__init__(self, parent, callback=callback, texts=texts,
                          objects=objects, categories=categories, *args, **kw)

if (__name__ == '__main__'):

  import Tkinter

  def callback(selected_indices, entry):

    print 'callback', selected_indices, entry

  root = Tkinter.Tk()

  menu = KeysymList(root, callback=callback)
  menu.grid()

  def getSelected():
    selected = menu.getSelected()
    print 'getSelected', selected

  def getSelectedInd():
    selected = menu.getSelectedIndex()
    print 'getSelectedInd', selected

  button = Tkinter.Button(root, text='getSelected', command=getSelected)
  button.grid()

  button = Tkinter.Button(root, text='getSelectedInd', command=getSelectedInd)
  button.grid()

  root.mainloop()
