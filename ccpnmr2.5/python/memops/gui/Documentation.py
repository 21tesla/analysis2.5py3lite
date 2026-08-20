
"""
======================COPYRIGHT/LICENSE START==========================

Documentation.py: <write function here>

Copyright (C) 2009 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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

class DocumentationException(Exception):
  pass

class WidgetDoc:

  def __init__(self, name, documentation=''):

    self.name = name
    self.documentation = documentation
    self.parentDoc = None

class ButtonDoc(WidgetDoc):

  pass

class ColumnDoc(WidgetDoc):

  def __init__(self, name, documentation='', isEditable=True):

    WidgetDoc.__init__(self, name, documentation)
    self.isEditable = isEditable

class TableDoc(WidgetDoc):

  def __init__(self, name, documentation='', columns=None):

    if columns is None:
      columns = ()

    WidgetDoc.__init__(self, name, documentation)
    self.columnDocs = []
    self.columnDict = {}

    for column in columns:
      column = tuple(column)
      self.addColumn(*column)

  def addColumn(self, name, documentation='', isEditable=True):

    columnDoc = ColumnDoc(name, documentation, isEditable)

    keys = set(self.columnDict.keys())
    n = 0
    while (name, n) in keys:
      n += 1
    self.columnDocs.append(columnDoc)
    self.columnDict[(name, n)] = columnDoc
    columnDoc.parentDoc = self

  def get(self, key):

    return self.columnDict.get(key)

  def __getitem__(self, key):

    return self.columnDict[key]

class ContainerDoc(WidgetDoc):

  def __init__(self, name=None, documentation='', widgetDocs=None):
    
    if widgetDocs is None:
      widgetDocs = ()

    WidgetDoc.__init__(self, name, documentation)

    self.parentDoc = None

    self.widgetDocs = []
    self.widgetDict = {}

    for widgetDoc in widgetDocs:
      self.addWidgetDoc(widgetDoc)

  def addWidgetDoc(self, widgetDoc):

    if not isinstance(widgetDoc, WidgetDoc):
      raise DocumentationException('widgetDoc = %s, must be a WidgetDoc')

    name = widgetDoc.name
    if name in self.widgetDict:
      raise DocumentationException('repeated widgetDoc name "%s"' % name)

    self.widgetDocs.append(widgetDoc)
    self.widgetDict[name] = widgetDoc
    widgetDoc.parentDoc = self

  def get(self, key):

    return self.widgetDict.get(key)

  def __getitem__(self, key):

    return self.widgetDict[key]

class TabDoc(ContainerDoc):

  pass

class PopupDoc(ContainerDoc):

  pass

cloneDoc = ButtonDoc('Clone', 'Clone popup window')
helpDoc = ButtonDoc('Help', 'Show popup help document')
closeDoc = ButtonDoc('Close', 'Close popup')
cancelDoc = ButtonDoc('Cancel', 'Cancel operation and close popup')

if __name__ == '__main__':

  popupDoc = PopupDoc(name='Edit Spectra',
    documentation='Allows the user to edit various details and parameters of the spectra within the current project',
    widgetDocs=(
      ButtonDoc('Delete Spectrum', 'Deletes the selected spectrum or spectra'),
      cloneDoc, helpDoc, closeDoc,
      TabDoc(name='Spectra',
        documentation='',
        widgetDocs=(
          TableDoc(
            name=None,
            columns=(
              ('#', 'Row number', False),
              ('Experiment', 'Experiment name', True),
              ('Spectrum', 'Spectrum name', True),
              ('Dims', 'Number of dimensions', False),
              ('Active Peaklist', 'The active peaklist, for example the peaklist into which picked peaks are put', False),
            ),
          ),
        ),
      ),
      TabDoc(name='Display Options',
      ),
    ),
  )
