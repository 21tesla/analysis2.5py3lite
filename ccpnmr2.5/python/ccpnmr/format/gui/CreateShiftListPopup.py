
"""
======================COPYRIGHT/LICENSE START==========================

CreateShiftListPopup.py: GUI window to set parameters for creating shift list from peak lists

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

from memops.universal.Io import joinPath

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.ScrolledListbox import ScrolledListbox
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.DataEntry import askString 
from memops.gui.MessageReporter import showInfo

from memops.universal.Util import returnFloat

from ccpnmr.format.general.Io import getHelpUrlDir

#
# Get script to create chemical shift list from peaklist(s)
#

from ccpnmr.format.process.makeChemShiftsFromPeakLists import createChemShifts

class CreateShiftListPopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'CreateShiftList.html')
  
  def __init__(self, parent, project):

    self.project = project
    self.peakLists = []
     
    BasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Create shift list from peak lists', modal=False, transient=True)

  def body(self, master):
  
    #
    # Peaklist setup
    #
    
    self.peakListDict = {}
    
    for experiment in self.project.currentNmrProject.sortedExperiments():
      for dataSource in experiment.sortedDataSources():
        for peakList in dataSource.sortedPeakLists():
          
          peakListLabel = "%s:%s:%d:%s" % (experiment.name,dataSource.name,peakList.serial,peakList.name)
          
          self.peakListDict[peakListLabel] = peakList
    
    peakListLabels = self.peakListDict.keys()
    peakListLabels.sort()
    
    #    
    # chemical shift list setup
    #
    
    self.shiftListDict = {}
    self.shiftListDict['None'] = None
    
    for shiftList in self.project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList'):
      
      shiftListLabel = "%d:%s" % (shiftList.serial,shiftList.name)
          
      self.shiftListDict[shiftListLabel] = shiftList
      
    shiftListLabels = self.shiftListDict.keys()
    shiftListLabels.sort()
      
    row = 0
    
    label = Label(master, text= "Generation of chemical shift list")
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.W)
    
    row += 1
    
    label = Label(master, text= "Peak lists:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.peakListBox = ScrolledListbox(master, width = 50, height = 5, selectmode = Tkinter.MULTIPLE,
                                initial_list = peakListLabels)
    self.peakListBox.grid(row=row, column=1, sticky=Tkinter.EW)

    row += 1
    
    label = Label(master, text= "Use existing shift list:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.shiftListSelect = PulldownMenu(master, entries = shiftListLabels)
    self.shiftListSelect.grid(row=row, column=1, sticky=Tkinter.EW)

    row += 1
    
    label = Label(master, text= "Use multiple assignments:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.useAllContribs = CheckButton(master)
    self.useAllContribs.grid(row=row, column=1, sticky=Tkinter.W)
    
    #
    # Setup the default shift error per nucleus
    #
    
    self.defaultShiftError = {}

    for (nucl,text,defValue) in [('1H','proton','0.002'),('13C','carbon','0.1'),('15N','nitrogen','0.1')]:

      row += 1

      label = Label(master, text= "Default %s shift error:" % text)
      label.grid(row=row, column=0, sticky=Tkinter.W)

      self.defaultShiftError[nucl] = Entry(master, text = defValue)
      self.defaultShiftError[nucl].grid(row=row, column=1, sticky=Tkinter.W)

    row += 1

    texts = [ 'Create shift list' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, columnspan = 2, column=0)

  def apply(self):
    
    selectedPeakLists = self.peakListBox.getSelectedItems()
    
    for peakListLabel in selectedPeakLists:
    
      self.peakLists.append(self.peakListDict[peakListLabel])
    

    chemShiftList = self.shiftListDict[self.shiftListSelect.getSelected()]
    
    if not chemShiftList:
    
      chemShiftListName = askString("Enter chemical shift list name","New chemical shift list name",'',self)
      
    else:
      
      chemShiftListName = chemShiftList.name
    
    useAllContribs = self.useAllContribs.isSelected()
    
    defaultShiftError = {}
    for nucl in self.defaultShiftError.keys():
      defaultShiftError[nucl] = returnFloat(self.defaultShiftError[nucl].get())
    
    listCreated = createChemShifts(self.peakLists, guiParent = self.parent, multiDialog = self.parent.multiDialog, shiftList = chemShiftList, useAllContribs = useAllContribs, defaultShiftError = defaultShiftError, shiftListName = chemShiftListName)
  
    if listCreated:
      showInfo("Success","Succesfully created a chemical shift list from peaklist(s)")
          
    return True

