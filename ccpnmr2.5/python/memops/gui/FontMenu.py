"""
======================COPYRIGHT/LICENSE START==========================

FontMenu.py: <write function here>

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

# import Tkinter as tk
from memops.gui.Menu import Menu

fontNames = ('Courier','Helvetica','Lucida','System','Times')
BOLD = 'bold'
ITALIC = 'italic'


class FontMenu(Menu):

    def __init__(self, parent, setFunc, sizes=(8, 10, 12, 14), doBold=1, doItalic=1, doBoldItalic=0, *args, **kw):

        Menu.__init__(self, parent, *args, **kw)

        # NOTE:ED potentially read all the fonts available on the system
        #       restricted to one-word names as must be stored as a string (tuple as nicer format)
        # fontNames = tkFont.families()
        # filter those that contain spaces
        # fontNames = [fnt for fnt in fontNames if len(fnt.split()) == 1]

        subMenus = {}
        for fontName in fontNames:
            subMenus[fontName] = Menu(self, tearoff=0)

            for s in sizes:
                spec = '%s %d' % (fontName, s)
                label = '%dpt' % s
                # fontSpec = (fontName, s)
                subMenus[fontName].add_command(label=label, font=spec, command=lambda fnt=spec: setFunc(fnt))

            if doBold:
                for s in sizes:
                    spec = '%s %d %s' % (fontName, s, BOLD)
                    label = '%dpt %s' % (s, BOLD)
                    # fontSpec = (fontName, s, BOLD)
                    subMenus[fontName].add_command(label=label, font=spec, command=lambda fnt=spec: setFunc(fnt))

            if doItalic:
                for s in sizes:
                    spec = '%s %d %s' % (fontName, s, ITALIC)
                    label = '%dpt %s' % (s, ITALIC)
                    # fontSpec = (fontName, s, ITALIC)
                    subMenus[fontName].add_command(label=label, font=spec, command=lambda fnt=spec: setFunc(fnt))

            if doBoldItalic:
                for s in sizes:
                    spec = '%s %d %s %s' % (fontName, s, BOLD, ITALIC)
                    label = '%dpt %s %s' % (s, BOLD, ITALIC)
                    # fontSpec = (fontName, s, BOLD, ITALIC)
                    subMenus[fontName].add_command(label=label, font=spec, command=lambda fnt=spec: setFunc(fnt))

            self.add_cascade(label=fontName, shortcut=fontName[0], menu=subMenus[fontName])
