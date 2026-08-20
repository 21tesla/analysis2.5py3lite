
"""
======================COPYRIGHT/LICENSE START==========================

WriteBmrbChemShiftDepPopup.py: GUI window for writing BMRB chemical shift deposition file

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

from memops.universal.Io import joinPath

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.MessageReporter import showInfo, showError

from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup

from ccpnmr.format.general.Io import getHelpUrlDir

#
# Get script to write BMRB deposition file
#

from ccpnmr.format.process.writeBmrbChemShiftDeposition import writeBmrbChemShiftDeposition

class WriteBmrbChemShiftDepPopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'WriteBmrbChemShiftDep.html')

  def __init__(self, parent, project):

    self.project = project
    self.guiParent = parent
    self.selectedFormats = []
    self.defaultText = 'Select file'
        
    self.chainDict = {}
    self.chainList = []
    
    self.shiftListDict = {}
    self.shiftLists = []
    
    for shiftList in self.project.currentNmrProject.findAllMeasurementLists(className = 'ShiftList'):
      if shiftList.measurements:
        label = str(shiftList.serial) + ':' + str(shiftList.name)
        self.shiftListDict[label] = shiftList
        self.shiftLists.append(label)
        
    if not self.shiftLists:
    
      showError('Error','No shift lists available!')
      return
    
    BasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Write BMRB chemical shift deposition file', modal=False, transient=True)

  def body(self, master):
        
    row = 0
    
    label = Label(master, text= "BMRB chemical shift deposition file writer.")
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.W)

    row += 1

    label = Label(master, text= "Shift lists:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.shiftListSelect = PulldownMenu(master, entries = self.shiftLists, callback = self.setChainList, do_initial_callback = False)
    self.shiftListSelect.grid(row=row, column=1, sticky=Tkinter.EW)
    
    row += 1
    
    label = Label(master, text= "Chains (only one per file):")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.chainListSelect = PulldownMenu(master, entries = self.chainList)
    self.chainListSelect.grid(row=row, column=1, sticky=Tkinter.EW)
    self.setChainList(0,self.shiftLists[0])
    
    row += 1
    
    label = Label(master, text= "Chemical shift deposition file:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.fileButton = Tkinter.Button(master, text = self.defaultText, command = self.selectFile)
    self.fileButton.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1

    texts = [ 'Write file' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, columnspan = 2, column=0)

  def setChainList(self,shiftListIndex,shiftListLabel):
  
    #
    # Has to be done... very slow this.
    #
  
    shiftList = self.shiftListDict[shiftListLabel]
     
    self.chainDict = {}
    self.chainList = []

    for shift in shiftList.measurements:
      
      if shift.resonance:

        rs = shift.resonance.resonanceSet

        if rs:
          atom = rs.findFirstAtomSet().findFirstAtom()
          chain = atom.residue.chain

          if chain not in self.chainDict.values():
            label = chain.molSystem.code + ":'" + chain.code + "'"
            self.chainDict[label] = chain
            self.chainList.append(label)
                   
    self.chainListSelect.clearMenuItems()
    
    if self.chainList:
      self.chainListSelect.setup(self.chainList,0)

  def selectFile(self):
   
    fileName = self.fileButton.__getitem__('text')
    
    if fileName == self.defaultText:
      fileName = 'bmrb.csdep'
    
    popup = FormatFilePopup(self, file = fileName, component = 'csdep', format = 'bmrb')
    
    if popup.fileSelected:
    
      self.fileButton.config(text = popup.file)      
      popup.destroy()

  def apply(self):
    
    shiftListLabel = self.shiftListSelect.getSelected()
    shiftList = self.shiftListDict[shiftListLabel]

    try:
      chainLabel = self.chainListSelect.getSelected()
      chain = self.chainDict[chainLabel]
    except:
      showError("No chains","No chains were present or selected. Try running linkResonances first.")
      return False
    
    fileName = self.fileButton.__getitem__('text')
    
    if fileName == self.defaultText:
      return False
    
    fileCreated = writeBmrbChemShiftDeposition(self.guiParent,chain,shiftList,fileName)
    
    if fileCreated:
      showInfo("Success","Succesfully wrote chemical shift deposition file")
    else:
      showError("Not written","Error writing file %s. File not written" % fileName)
      return False
          
    return True

