
"""
======================COPYRIGHT/LICENSE START==========================

AutoAssignProjectInfoPopup.py: GUI window getting info to write AutoAssign project

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnFloat

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.MessageReporter import showInfo, showError
from memops.gui.Separator import Separator

from ccpnmr.format.gui.ObjectEditPopup import createSelection

from ccpnmr.format.general.Io import getHelpUrlDir


class AutoAssignProjectInfoPopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'AutoAssignProjectInfo.html')

  def __init__(self, parent, project,autoAssignProjectFile):

    self.project = project
    self.guiParent = parent
    self.autoAssignProjectFile = autoAssignProjectFile
     
    self.validSelection = False
    
    self.peakLists = []
    
    for nmrExp in self.project.currentNmrProject.sortedExperiments():
      for ds in nmrExp.sortedDataSources():
        for pl in ds.sortedPeakLists():
          self.peakLists.append(pl)
    
    if not self.peakLists:
      raise('No peak lists available!')
    
    (self.peakNameList,self.peakNameDict) = createSelection(self.peakLists)
    self.peakNameList.insert(0,'None')
    self.peakNameDict['None'] = None
      
    BasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Write AutoAssign project file', modal=True, transient=True)

  def body(self, master):
    
    self.geometry('600x400')

    columnspan = 5
    
    row = 0
    
    label = Label(master, text= "AutoAssign project file setup window.")
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.W)
    
    #
    # Set up the spectra
    #
    
    row += 1
    
    label = Label(master, text= "Spectrum type")
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.W)

    label = Label(master, text= "CCPN peak list")
    label.grid(row=row, column=2, columnspan = 2, sticky=Tkinter.W)

    label = Label(master, text= "Phase") # TODO What's this?
    label.grid(row=row, column=4, sticky=Tkinter.W)
    
    self.spectrumSelectWidget = {}
    self.spectrumPhaseWidget = {}

    for i in range(len(self.autoAssignProjectFile.defaultSpectra)):

      row += 1
    
      specType = self.autoAssignProjectFile.defaultSpectra[i][0]
      
      if i == 0:
        labelText = "%s (root)" % specType
      else:
        labelText = "%s" % specType

      label = Label(master, text= labelText)
      label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.W)

      self.spectrumSelectWidget[specType] = PulldownMenu(master, entries = self.peakNameList)
      self.spectrumSelectWidget[specType].grid(row=row, column=2, columnspan = 2, sticky=Tkinter.W)
    
      self.spectrumPhaseWidget[specType] = Entry(master, text = "", width = 4)
      self.spectrumPhaseWidget[specType].grid(row=row, column=4, sticky=Tkinter.W)

    
    #
    # Make a break...
    #

    row += 1
    
    separator = Separator(master,height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
    
    #
    # Set up the tolerances
    #

    row += 1
    
    label = Label(master, text= "Atom tolerances")
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
    
    toleranceAtoms = self.autoAssignProjectFile.toleranceAtoms
    toleranceDefault = self.autoAssignProjectFile.toleranceDefault
    
    self.toleranceWidget = {}
    
    for i in range(0,len(toleranceAtoms),2):
      row += 1
      for j in range(2):
        curPos = i + j
        
        if curPos == len(toleranceAtoms):
          break
        
        toleranceAtom = toleranceAtoms[curPos]
        
        label = Label(master, text= toleranceAtom)
        label.grid(row=row, column=2 * j, sticky=Tkinter.E)
        
        self.toleranceWidget[toleranceAtom] = Entry(master, text = str(toleranceDefault[curPos]), width = 6)
        self.toleranceWidget[toleranceAtom].grid(row=row, column= 1 + 2 * j, sticky=Tkinter.W)

    #
    # Make a break...
    #

    row += 1
    
    separator = Separator(master,height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)

    row += 1
    
    label = Label(master, text= "Will save the peak files in the project file directory.", fg = 'red')
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
    
    row += 1

    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, columnspan = columnspan, column=0)

  def apply(self):
    
    rootSpecType = self.autoAssignProjectFile.defaultSpectra[0][0]
    rootPeakList = self.peakNameDict[self.spectrumSelectWidget[rootSpecType].getSelected()]
      
    if not rootPeakList:
      showError("No root spectrum","Can't continue: need a root spectrum...")
      return False
      
    #
    # Organize peak/spectrum info
    #
    
    self.peakListInfo = []
            
    for i in range(0,len(self.autoAssignProjectFile.defaultSpectra)):
    
      specType = self.autoAssignProjectFile.defaultSpectra[i][0]
      peakList = self.peakNameDict[self.spectrumSelectWidget[specType].getSelected()]
      phase = self.spectrumPhaseWidget[specType].get()
      
      if peakList:
        self.peakListInfo.append((specType,peakList,phase))

    #
    # Organize tolerance info
    #
    
    self.tolerances = {}
    
    toleranceAtoms = self.autoAssignProjectFile.toleranceAtoms
    toleranceDefault = self.autoAssignProjectFile.toleranceDefault
    
    for i in range(len(toleranceAtoms)):
        
      toleranceAtom = toleranceAtoms[i]
      value = returnFloat(self.toleranceWidget[toleranceAtom].get())
      
      if not value:
        value = toleranceDefault[i]
      
      self.tolerances[toleranceAtom] = value

    self.validSelection = True
    
    return True
