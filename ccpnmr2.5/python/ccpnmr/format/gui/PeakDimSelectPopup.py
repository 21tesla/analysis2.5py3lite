
"""
======================COPYRIGHT/LICENSE START==========================

PeakDimSelectPopup.py: GUI window for matching experiment dimensions in data model to exported peak list dimensions

Copyright (C) 2005-2011 Wim Vranken (Vrije Universiteit Brussel)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
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

- contact Wim Vranken (wim@ebi.ac.uk)
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
import Tkinter

from memops.universal.Io import joinPath

from memops.universal.Util import returnInt
from ccpnmr.format.general.Io import getHelpUrlDir

from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Label import Label
from memops.gui.Util import createHelpButtonList

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class PeakDimSelectPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'PeakDimSelect.html')
   
  def __init__(self, parent, peakList, format, dataDimRefs = None, order = None):
  
    self.peakList = peakList
    self.dataSource = self.peakList.dataSource
    self.format = format
    
    if dataDimRefs:
      self.dataDimRefs = dataDimRefs

    else:
      self.dataDimRefs = []

    numValidDims = 0
    for dataDim in self.peakList.dataSource.dataDims:
      if dataDim.className != 'SampledDataDim':
        numValidDims += 1
    
    self.numDim = numValidDims
    self.order = order[self.numDim]
    
    TemporaryBasePopup.__init__(self,parent = parent, title = "Project '%s': " % self.peakList.root.name + "DataDimRef selection for PeakDims", modal = False, transient=True)
 
  def body(self, master):

    #
    # Initialize
    #
    
    #
    # Get relevant dataDimRefs AND chemical shift ranges for each
    #

    self.dataDimRefDict = {}
    dataDimRefList = []
    dataDimRefSelection = self.numDim * ['']
    
    for peak in self.peakList.sortedPeaks():
  
      for peakDim in peak.sortedPeakDims():
      
        dataDim = peakDim.dataDim
        
        if dataDim.className == 'SampledDataDim':
          continue
        
        dataDimRef = peakDim.dataDimRef
        
        if dataDim.dim == 1:
          addInfo = " (acqu)"
        else:
          addInfo = ""
            
        isotopeString = '/'.join(dataDimRef.expDimRef.isotopeCodes)
        
        selectionString = "Dim %d, nucl %s%s" % (dataDim.dim,isotopeString,addInfo)

        if dataDimRefList.count(selectionString) == 0:
          dataDimRefList.append(selectionString)
          self.dataDimRefDict[selectionString] = dataDimRef
          
          if dataDimRef in self.dataDimRefs:
            dataDimRefSelection[self.dataDimRefs.index(dataDimRef)] = selectionString
          else:
            dataDimRefSelection[(peakDim.dim - 1)] = selectionString
    
    #
    # Popup info
    #
    
    #
    # Header labels
    #

    row = 0
    
    label = Label(master, text= "%s peak dim" % self.format)
    label.grid(row=row, column=0, sticky=Tkinter.EW)
    
    label = Label(master, text= "PeakDim selection")
    label.grid(row=row, column=1, sticky=Tkinter.EW)

    #
    # Selection per peakDim
    #
    
    self.dataDimRefMenu = []
    
    for peakDim in range(0,self.numDim):

      row = row + 1

      label = Label(master, text= str(peakDim))
      label.grid(row=row, column=0, sticky=Tkinter.EW)
      
      peakDimIndex = self.order[peakDim]

      self.dataDimRefMenu.append(PulldownMenu(master, entries = dataDimRefList, selected_index = dataDimRefList.index(dataDimRefSelection[peakDimIndex])))
      self.dataDimRefMenu[-1].grid(row=row, column=1, sticky=Tkinter.W, ipadx = 20)
  
    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 2)
   

  def apply(self):
    
    self.dataDimRefs = []
    
    for peakDim in range(0,self.numDim):
      
      dataDimRef = self.dataDimRefDict[self.dataDimRefMenu[peakDim].getSelected()]
      
      if self.dataDimRefs.count(dataDimRef) > 0:
        self.dataDimRefs = []
        return False
      
      else:
        self.dataDimRefs.append(dataDimRef)
    
    return True
