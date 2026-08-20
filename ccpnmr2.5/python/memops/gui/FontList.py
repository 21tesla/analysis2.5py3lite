
"""
======================COPYRIGHT/LICENSE START==========================

FontList.py: Part of the CcpNmr software

Copyright (C) 2010 Wayne Boucher and Tim Stevens (University of Cambridge)

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
from memops.gui.Frame import Frame
from memops.gui.PulldownList import PulldownList

from memops.universal import Output

tkNames = [ 'Helvetica', 'Times', 'Courier' ]
tkSizes = [ 6, 8, 10, 12, 14, 16, 18, 20, 22, 24 ]

glNames = tkNames + ['Roman', 'MonoRoman']
glSizes = {
  'Helvetica': [ 10, 12, 18 ],
  'Times': [ 10, 24 ],
  'Courier': [ 13, 15 ],
  'Roman': tkSizes,
  'MonoRoman': tkSizes,
}

printNames = Output.font_names
printSizes = Output.font_sizes

fontModes = ['OpenGL', 'Tk', 'Print']

class FontList(PulldownList):
 
  def __init__(self, parent, selected=None,
               mode='Tk', callback=None, extraTexts=None, *args, **kw):
 
    assert mode in fontModes, 'mode = %s' % mode

    self.mode = mode
    
    if not selected:
      selected = 'Helvetica 10'

    if not extraTexts:
      extraTexts = []

    PulldownList.__init__(self, parent, callback, *args, **kw)
 
    texts = extraTexts
    cats = len(texts)*[None]
 
    if mode == 'OpenGL':
      for name in glNames:
        for size in glSizes[name]:
          text = '%s %d' % (name, size)
          texts.append(text)
          cats.append(name)
      
    elif mode == 'Tk':
      for name in tkNames:
        for size in tkSizes:
          text = '%s %d' % (name, size)
          texts.append(text)
          cats.append(name)

    else: # mode == 'Print':
      for name in printNames:
        for size in printSizes:
          text = '%s %d' % (name, size)
          texts.append(text)
          cats.append(name)

    if selected in texts:
      index = texts.index(selected)
    else:
      index = 0

    PulldownList.setup(self, texts, texts, index, categories=cats)

if __name__ == '__main__':

  def myCallback(font):
    print 'myCallback', font

  import Tkinter
  r = Tkinter.Tk()
  f = FontList(r, isOpenGL=True, callback=myCallback)
  f.grid()
  r.mainloop()
