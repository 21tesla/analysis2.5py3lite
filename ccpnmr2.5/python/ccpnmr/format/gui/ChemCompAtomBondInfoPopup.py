
"""
======================COPYRIGHT/LICENSE START==========================

ChemCompAtomBondInfoPopup.py: GUI window for defining chemComp atom and bond info

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

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.CheckButton import CheckButton
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Util import createHelpButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.Separator import Separator
from memops.gui.MessageReporter import showWarning

from ccp.general.Constants import standardBackboneAtoms

from ccpnmr.format.general.Io import getHelpUrlDir

class ChemCompAtomBondInfoPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ChemCompAtomBondInfo.html')

  def __init__(self, parent, chemCompInfo, chemCompAtoms, chemCompBonds):
   
    self.chemCompInfo = chemCompInfo
    self.chemCompAtoms = chemCompAtoms
    self.chemCompBonds = chemCompBonds
    
    self.obligatoryAtoms = []
    
    for molTypeTuple in standardBackboneAtoms.keys():
      if self.chemCompInfo[0][1] in molTypeTuple:
        for atomName in standardBackboneAtoms[molTypeTuple]:
          if atomName not in self.obligatoryAtoms:
            self.obligatoryAtoms.append(atomName)
      
    
    self.nonEntryAttributes = {
    
      'chirality':   (PulldownMenu,['R','S','unknown']),
      'waterExchangeable': (CheckButton,False)
    
    }
        
    TemporaryBasePopup.__init__(self, parent=parent, title='ChemComp atom and bond creation', modal=False, transient=True)

  def body(self, master):
      
    #
    # Popup window
    #
    
    self.columnspan = len(self.chemCompAtoms[0])

    row = 0
    label = Label(master, text= "Molecule type '%s', ccp code '%s'" % (self.chemCompInfo[0][1],self.chemCompInfo[1][1]))
    label.grid(row=row, column=0, columnspan = self.columnspan, sticky=Tkinter.EW)
    
    #
    # Warning message if necessary!
    #
    
    if self.obligatoryAtoms:
      row += 1
      
      label = Label(master, text= "Warning: have to define atoms: %s!" % string.join(self.obligatoryAtoms,', '), fg = 'red')
      label.grid(row=row, column=0, columnspan = self.columnspan, sticky=Tkinter.EW)
    
    #
    # Show top information...
    #
    
    row = row + 1
    colNum = 0
    self.frameWidth = 0
    
    for chemCompAtomInfo in self.chemCompAtoms[0]:
      
      attrName = chemCompAtomInfo[0]

      label = Label(master, text = attrName, width = len(attrName) + 2)
      label.grid(row=row, column=colNum, sticky=Tkinter.EW)
      
      self.frameWidth += len(attrName) + 2
      colNum += 1
    
    #
    # Separator
    #
    
    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = self.columnspan, sticky = Tkinter.EW)
      
    #
    # Separate frame for info (have to be able to scroll)
    #
    
    row += 1
    
    # THIS BIT TELLS MASTER TO CONFIGURE WINDOW ON INSIDE WIDGET!!
    master.grid_rowconfigure(row,weight = 1)
    for i in range(self.columnspan):
      master.grid_columnconfigure(i,weight = 1)
    
    self.atomFrame = ScrolledFrame(master, width = self.frameWidth, height = 300, doExtraConfig = False)
    self.atomFrame.grid(row=row, column=0, columnspan = self.columnspan, sticky=Tkinter.NSEW)
    self.atomFrameRow = row
    self.atomFrameMaster = master

    #
    # Separator
    #
    
    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = self.columnspan, sticky = Tkinter.EW)
          
    #
    # End bit...
    #
    
    row = row + 1
    texts = [ 'OK', 'Update' ]
    commands = [ self.ok, self.update ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = self.columnspan)
    
    self.setupAtomFrame()

  def setupAtomFrame(self,resetFrame = False, resetStatus = True):
    
    frameRow = 0
    self.atomWidgets = []
    
    if resetFrame:
      self.atomFrame = ScrolledFrame(self.atomFrameMaster, width = self.frameWidth, height = 300, doExtraConfig = False)
      self.atomFrame.grid(row=self.atomFrameRow, column=0, columnspan = self.columnspan, sticky=Tkinter.NSEW)    
    
    frame = self.atomFrame.frame
      
    #
    # Make fake info at top for col length...
    #
    
    colNum = 0
    
    for chemCompAtomInfo in self.chemCompAtoms[0]:
      
      attrName = chemCompAtomInfo[0]

      label = Label(frame, text = '', width = len(attrName) + 2)
      label.grid(row=frameRow, column=colNum, sticky=Tkinter.EW)
      
      colNum += 1
            
    #
    # Now atom info...
    #
    
    for chemCompAtom in self.chemCompAtoms:
      
      frameRow += 1
      colNum = 0
      
      self.atomWidgets.append([])

      for chemCompAtomInfo in chemCompAtom:

        attrName = chemCompAtomInfo[0]

        if attrName == 'bondedAtoms':

          bondedAtomText = ""
          
          for otherAtom in chemCompAtomInfo[1]:
            bondedAtomText += otherAtom[0][1] + ','
          
          label = Label(frame, text = bondedAtomText[:-1])
          label.grid(row=frameRow, column=colNum, sticky=Tkinter.EW)
            
        else:

          if attrName in self.nonEntryAttributes:

            widgetInfo = self.nonEntryAttributes[attrName]

            if widgetInfo[0] == PulldownMenu:

              self.atomWidgets[-1].append(PulldownMenu(frame, entries = widgetInfo[1], selected_index = widgetInfo[1].index(chemCompAtomInfo[1])))

            elif widgetInfo[0] == CheckButton:

              self.atomWidgets[-1].append(CheckButton(frame, selected = chemCompAtomInfo[1]))

          else:

            text = chemCompAtomInfo[1]
            if not text:
              text = ''

            self.atomWidgets[-1].append(Entry(frame, text = text, width = 5))

          self.atomWidgets[-1][-1].grid(row=frameRow, column=colNum)     
      
        colNum += 1

  def update(self):
  
    self.apply(checkAtoms = False)
    self.setupAtomFrame()

  def apply(self, checkAtoms = True):
    
    atomNames = []
    duplicateNames = []
    
    if checkAtoms:
      obligatoryAtoms = self.obligatoryAtoms[:]
  
    for atomNum in range(0,len(self.chemCompAtoms)):
    
      chemCompAtom = self.chemCompAtoms[atomNum]
      widgets = self.atomWidgets[atomNum]
      
      for widgetNum in range(0,len(widgets)):

        attrName = chemCompAtom[widgetNum][0]
        
        if attrName in self.nonEntryAttributes:

          widgetInfo = self.nonEntryAttributes[attrName]

          if widgetInfo[0] in [PulldownMenu,CheckButton]:

            value = widgets[widgetNum].getSelected()
            
            if widgetInfo[0] == CheckButton:
              if value:
                value = True
              else:
                value = False

        else:

          value = widgets[widgetNum].get()
          
          if value == '':
            value = None

        chemCompAtom[widgetNum] = (attrName,value)

        if attrName == 'name':
          
          if value in atomNames:
            if value not in duplicateNames:
              duplicateNames.append(value)
          else:
            atomNames.append(value)
          
          if checkAtoms:
            if value in obligatoryAtoms:
              obligatoryAtoms.pop(obligatoryAtoms.index(value))
    
    if duplicateNames:
      showWarning(self,"Duplicate atom names: %s - names have to be unique!" % string.join(duplicateNames,', '))
      return False
    
    if checkAtoms and obligatoryAtoms:
      showWarning(self,"Obligatory atoms: %s still missing - please rename accordingly." % string.join(obligatoryAtoms,', '))
      return False
      
    return True
  
