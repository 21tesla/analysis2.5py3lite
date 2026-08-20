
"""
======================COPYRIGHT/LICENSE START==========================

DataDimRefSelectByPeakDimPopup.py: GUI window for matching original peak list dimensions to experiment dimensions in data model

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

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

import string


from ccpnmr.format.general.Io import getHelpUrlDir
from memops.universal.Constants import dirsep
from memops.universal.Io import joinPath

from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Label import Label
from memops.gui.BasePopup import BasePopup
from memops.gui.Util import createHelpButtonList

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class DataDimRefSelectByPeakDimPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'DataDimRefSelectByPeakDim.html')

  def __init__(self, parent, dataSource, peaks, order = None):
  
    self.dataSource = dataSource
    self.numDim = len(self.dataSource.dataDims)
    self.peaks = peaks
    self.dataDimRefs = []
    self.order = order[self.numDim]
    
    TemporaryBasePopup.__init__(self,parent = parent, title = "Project '%s': " % dataSource.root.name + "DataDimRef selection for PeakDims", modal = False, transient=True)
 
  def body(self, master):

    #
    # Initialize
    #
    
    #
    # Get chem shift range
    #
    
    peakChemShiftRange = []
    
    if not self.peaks:
      return
    elif hasattr(self.peaks[0],'ppm'):
      peakUnit = 'ppm'
    elif hasattr(self.peaks[0],'point'):
      peakUnit = 'point'
    else:
      return
    
    for peak in self.peaks:
      
      peakInfo = getattr(peak,peakUnit)

      for dim in range(0,len(peakInfo)):
      
        if len(peakChemShiftRange) < dim + 1:
         peakChemShiftRange.append([999999,-999999])
    
        if peakChemShiftRange[dim][0] > peakInfo[dim]:
          peakChemShiftRange[dim][0] = peakInfo[dim]

        if peakChemShiftRange[dim][1] < peakInfo[dim]:
          peakChemShiftRange[dim][1] = peakInfo[dim]
         
    #
    # Make datadimref lists
    #
    
    dataDimRefList = []
    self.dataDimRefDict = {}
    dataDimRefSelection = []
    
    for i in self.order:
      
      dataDim = self.dataSource.sortedDataDims()[i]
      
      if dataDim.dim == 1:
        addInfo = " (acqu)"
      else:
        addInfo = ""
    
      for dataDimRef in dataDim.dataDimRefs:
        
        # Note: key for dataDimRef is expDimRef!
        
        isotopeString = string.join(dataDimRef.expDimRef.isotopeCodes,dirsep)
        
        selectionString = "Dim %d, nucl %s, refValue %s%s" % (dataDim.dim,isotopeString,dataDimRef.refValue,addInfo)
        
        # In principle this could go wrong if one expDim has two expDimRefs with same isotopelist
        # (but what would be the use? Is this a good constraint?)
        dataDimRefList.append(selectionString)
        self.dataDimRefDict[selectionString] = dataDimRef
        
      dataDimRefSelection.append(dataDimRefList[-1])
    
    #
    # Add on a None dimension if necessary
    #
    
    if self.numDim < len(peakChemShiftRange):
      dataDimRefSelection.append('None')
      dataDimRefList.append('None')
      self.dataDimRefDict['None'] = None
    
    #
    # Popup info
    #
    
    #
    # Header labels
    #

    row = 0
    
    label = Label(master, text= "Peak dim")
    label.grid(row=row, column=0, sticky=Tkinter.EW)
    
    label = Label(master, text= "Chem. shift range")
    label.grid(row=row, column=1, sticky=Tkinter.EW)
    
    # TODO: add in peak dim column headers if known (for nmrView for example)
    
    label = Label(master, text= "DataDimRef selection")
    label.grid(row=row, column=2, sticky=Tkinter.EW)

    #
    # Selection per peakDim
    #
    
    self.dataDimRefMenu = []
    
    for peakDim in range(0,len(peakChemShiftRange)):

      row = row + 1

      label = Label(master, text= str(peakDim))
      label.grid(row=row, column=0, sticky=Tkinter.EW)

      label = Label(master, text= "%.3f - %.3f %s" % (peakChemShiftRange[peakDim][0], peakChemShiftRange[peakDim][1], peakUnit))
      label.grid(row=row, column=1, sticky=Tkinter.EW)
      
      if peakDim < len(dataDimRefList):
        selectedIndex = dataDimRefList.index(dataDimRefSelection[peakDim])
      else:
        selectedIndex = len(dataDimRefList) - 1

      self.dataDimRefMenu.append(PulldownMenu(master, entries = dataDimRefList, selected_index = selectedIndex))
      self.dataDimRefMenu[-1].grid(row=row, column=2, sticky=Tkinter.E, ipadx = 20)
 
 
    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 3)
   

  def apply(self):
    
    for peakDim in range(0,len(self.dataDimRefMenu)):
      
      dataDimRef = self.dataDimRefDict[self.dataDimRefMenu[peakDim].getSelected()]
      
      if dataDimRef and self.dataDimRefs.count(dataDimRef) > 0:
        self.dataDimRefs = []
        return False
      
      self.dataDimRefs.append(dataDimRef)
    
    if (len(self.dataDimRefs) - self.dataDimRefs.count(None)) != self.numDim:
      self.dataDimRefs = []
      return False
    
    return True
