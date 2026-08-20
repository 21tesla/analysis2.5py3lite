
"""
======================COPYRIGHT/LICENSE START==========================

SequenceChemCompSelect.py: GUI window for setting the chemComp/chemCompVar for a sequence

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

from memops.universal.Util import returnStrings
from memops.universal.Util import returnInt

from memops.gui.Label import Label
from memops.gui.CheckButton import CheckButton
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.Separator import Separator
from memops.gui.MessageReporter import showError
from memops.gui.BasePopup import BasePopup

#from ccp.api.nmr import Nmr, NmrConstraint

from ccpnmr.format.gui.ChemCompSelectionPopup import ChemCompSelectionPopup

from ccpnmr.format.general.Io import getHelpUrlDir

class SequenceChemCompSelect(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'SequenceChemCompSelect.html')

  def __init__(self, parent, project, molName, createMoleculeDict, localCreateMoleculeDict):

    self.molName = molName
    self.project = project
    self.createMoleculeDict = createMoleculeDict
    self.localCreateMoleculeDict = localCreateMoleculeDict
    
    self.hasChanged = False
    
    self.molecule = self.localCreateMoleculeDict[self.molName][0]
    
    self.deleteText = 'Delete'
    self.okText = 'OK'
    self.startText = 'Insert at start'
    self.endText = 'Append to end'
    
    self.newMoleculeList = []
    
    for molName in self.localCreateMoleculeDict.keys():
      if molName != self.molName and not self.localCreateMoleculeDict[molName][0]:
        self.newMoleculeList.append(molName)
       
    # modal = true means that it won't continue unless this one returns value
    BasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name +'Set chemComps for sequence', modal=True, transient=True)

  def body(self, master):

    self.geometry('700x500')

    #
    # Setup header 
    #
    
    row = 0

    label = Label(master, text = 'Original', width = 10)
    label.grid(row=row, column=0, sticky=Tkinter.EW)

    label = Label(master, text = 'Modify ChemCompVar', width = 35)
    label.grid(row=row, column=1, sticky=Tkinter.EW)
    
    label = Label(master, text = 'Status', width = 20)
    label.grid(row=row, column=2, sticky=Tkinter.EW)
    
    label = Label(master, text = '', width = 5)
    label.grid(row=row, column=3, sticky=Tkinter.EW)

    row += 1

    label = Label(master, text = 'code', width = 10)
    label.grid(row=row, column=0, sticky=Tkinter.EW)

    label = Label(master, text = 'ccpCode (molType): linking, descriptor', width = 35)
    label.grid(row=row, column=1, sticky=Tkinter.EW)
    
    label = Tkinter.Button(master, text = 'Add new molecule', command = self.addNewMolecule)
    label.grid(row=row, column=2)
    
    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = 3, sticky = Tkinter.EW)
   
    row += 1
    
    # THIS BIT TELLS MASTER TO CONFIGURE WINDOW ON INSIDE WIDGET!!
    master.grid_rowconfigure(row,weight = 1)
    for i in range(4):
      master.grid_columnconfigure(i,weight = 1)
    self.sequenceFrame = ScrolledFrame(master, width = 70, height = 300, doExtraConfig = False)
    self.sequenceFrame.grid(row=row, column=0, columnspan = 4, sticky=Tkinter.NSEW)
    self.sequenceFrameRow = row
    self.sequenceFrameMaster = master
    
    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = 3, sticky = Tkinter.EW)

    row += 1

    texts = [ 'Change', 'Change and quit']
    commands = [ self.updateMolDict, self.ok]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Cancel', help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)

    self.setupSequenceFrame()    

  def setupSequenceFrame(self, resetFrame = False, resetStatus = True):
    
    frameRow = 0
    x = y = 0
    
    if resetFrame:
      self.sequenceFrame = ScrolledFrame(self.sequenceFrameMaster, width = 70, height = 300, doExtraConfig = False)
      self.sequenceFrame.grid(row=self.sequenceFrameRow, column=0, columnspan = 4, sticky=Tkinter.NSEW)    
    
    frame = self.sequenceFrame.frame

    #
    # Just use the first one as reference!
    #

    sequences = self.localCreateMoleculeDict[self.molName][1]
    (self.origSequence,self.origSequenceList) = sequences[0]
    
    if len(sequences) > 1:
      multipleSequences = True
    else:
      multipleSequences = False
      
    #
    # Set info...
    #

    label = Label(frame, text = '', width = 10)
    label.grid(row=frameRow, column=0, sticky=Tkinter.EW)

    label = Label(frame, text = '', width = 35)
    label.grid(row=frameRow, column=1, sticky=Tkinter.EW)
    
    label = Label(frame, text = '', width = 20)
    label.grid(row=frameRow, column=2, sticky=Tkinter.EW)
    
    frameRow += 1
    
    self.ccvButtons = []
    
    if resetStatus or resetFrame:
      self.statusList = []
      self.moveStatus = {}

      #
      # Set up a list of status entries
      #

      statusEntries = [self.okText,self.deleteText]

      for newMolName in self.newMoleculeList:
        for appendCode in [self.startText,self.endText]:
          statusEntries.append('%s of %s' % (appendCode,newMolName))
          self.moveStatus[statusEntries[-1]] = (newMolName,appendCode)
      
    #
    # Start creating all the objects...
    #

    for seqIndex in range(len(self.origSequenceList)):

      (seqEl,chemComp,ccVar) = self.origSequenceList[seqIndex]
      
      #
      # Set the original sequence code info...
      #
      
      seqCode = seqEl.seqCode

      if hasattr(seqEl,'code3Letter') and seqEl.code3Letter:
        seqText = seqEl.code3Letter
      else:
        seqText = seqEl.code1Letter

      label = Label(frame, text = "%s-%d" % (seqText,seqCode))
      label.grid(row = frameRow, column = 0, sticky = Tkinter.EW)
      
      #
      # Set the chemCompVar info
      #
      
      keywds = {}

      if chemComp:
        textCode = "%s (%s)" % (chemComp.ccpCode,chemComp.molType)
        
        if ccVar:
          textCode += ": %s,%s" % (ccVar.linking,ccVar.descriptor)
        else:
          textCode = 'Select chemCompVar for %s' % textCode
          keywds['fg'] = 'red'

      else:
        textCode = 'Select chemCompVar'
        keywds['fg'] = 'red'

      self.ccvButtons.append(Tkinter.Button(frame, text = textCode, command = lambda seqIndex = seqIndex:self.selectChemCompVar(seqIndex), **keywds))
      self.ccvButtons[-1].grid(row = frameRow, column = 1, sticky = Tkinter.EW)

      #
      # Set the selector...
      #
      
      if resetStatus or resetFrame:
        self.statusList.append(PulldownMenu(frame, entries = statusEntries, label_color = 'red'))
        self.statusList[-1].grid(row=frameRow, column=2, sticky=Tkinter.EW)
      
      frameRow += 1
      
    return True
  
  def addNewMolecule(self):
   
    molFormat = 'newMolecule_%d'
    
    if not hasattr(self,'molNames'):
      self.molNames = self.localCreateMoleculeDict.keys()
    
    for i in range(0,10):
      newMolName = molFormat % (i + 1)
      if newMolName not in self.molNames and not self.project.findFirstMolecule(name = newMolName):
        break
    
    self.molNames.append(newMolName)
    
    for statusPullDown in self.statusList:
      for appendCode in [self.startText,self.endText]:
        entryText = '%s of %s' % (appendCode,newMolName)
        statusPullDown.append(entryText)
        self.moveStatus[entryText] = (newMolName,appendCode)

  def updateMolDict(self):
    
    hasChanged = False
    actualSeqIndex = 0
    seqLength = len(self.origSequenceList)
    
    origSequenceList = self.origSequenceList[:]
  
    for seqIndex in range(0,seqLength):

      (seqEl,chemComp,ccVar) = origSequenceList[seqIndex]
    
      newStatus = self.statusList[seqIndex].getSelected()
      
      if newStatus != self.okText:
      
        hasChanged = True
      
        if self.localCreateMoleculeDict == self.createMoleculeDict:
          self.copyCreateMoleculeDict()

        self.origSequenceList.pop(actualSeqIndex)

        if newStatus != self.deleteText:

          (moveToMoleculeName,appendCode) = self.moveStatus[newStatus]
          
          if not self.localCreateMoleculeDict.has_key(moveToMoleculeName):
            self.localCreateMoleculeDict[moveToMoleculeName] = [None,[(None,[])]]

          sequences = self.localCreateMoleculeDict[moveToMoleculeName][1]
          
          for (tempSequence,tempSequenceList) in sequences:
            if appendCode == self.startText:
              tempSequenceList.insert(0,(seqEl,chemComp,ccVar))
            elif appendCode == self.endText:
              tempSequenceList.append((seqEl,chemComp,ccVar))
      
      else:
        actualSeqIndex += 1
       
    #
    # Check if changed and do some resets if so...
    #
    
    if hasChanged:
      self.hasChanged = True
      
      #
      # Reset the linking if appropriate
      #
      
      for molName in self.localCreateMoleculeDict.keys():
        if not self.localCreateMoleculeDict[molName][0]:
                  
          sequences = self.localCreateMoleculeDict[molName][1]
          
          for (origSequence,origSequenceList) in sequences:
            molTypes = []
            for (seqEl,chemComp,ccVar) in origSequenceList:
              if chemComp:
                if chemComp.molType not in molTypes:
                  molTypes.append(chemComp.molType)

            if len(molTypes) == 1 and molTypes[0] != 'other':
              seqLength = len(origSequenceList)
              for seqIndex in range(0,seqLength):
                (seqEl,chemComp,ccVar) = origSequenceList[seqIndex]
                if ccVar:
                  if seqIndex == 0:
                    linking = 'start'
                  
                  elif seqIndex == seqLength - 1:
                    linking = 'end'
                  
                  else:
                    linking = 'middle'

                  if ccVar.linking != linking:
                    ccVar = chemComp.findFirstChemCompVar(descriptor = ccVar.descriptor, linking = linking)
                    if ccVar:
                      origSequenceList[seqIndex] = (seqEl,chemComp,ccVar)

      self.setupSequenceFrame(resetFrame = True)
        
  def selectChemCompVar(self,seqIndex):
  
    hasChanged = False
  
    (seqEl,chemComp,ccVar) = self.origSequenceList[seqIndex]

    keywds = {}
    if chemComp:
      keywds['molTypeEntries'] = [chemComp.molType]  # TODO HAS TO CHANGE!!
      keywds['selectedChemComps'] = [chemComp]
    if ccVar:
      keywds['selectLinking'] = ccVar.linking
      
    if hasattr(seqEl,'code3Letter') and seqEl.code3Letter:
      keywds['origCode'] =  seqEl.code3Letter
    else:
      keywds['origCode'] =  seqEl.code1Letter

    popup = ChemCompSelectionPopup(self,self.project,chemCompEntries = ['ChemCompVar'], **keywds)
    self.wait_window(popup)
    
    if popup.frame.getSelectedChemComp():

      newChemCompVar = popup.frame.getSelectedChemComp()
      
      if not ccVar or ccVar != newChemCompVar:
        if self.localCreateMoleculeDict == self.createMoleculeDict:
          self.copyCreateMoleculeDict()
          
        self.origSequenceList[seqIndex] = (seqEl,newChemCompVar.chemComp,newChemCompVar)
        hasChanged = True
    
    if hasChanged:
      self.setupSequenceFrame(resetStatus = False)
      self.hasChanged = hasChanged
        
  def copyCreateMoleculeDict(self):
  
    self.localCreateMoleculeDict = {}
    
    for molName in self.createMoleculeDict.keys():
      self.localCreateMoleculeDict[molName] = [self.createMoleculeDict[molName][0],[]]
      
      for (sequence,sequenceList) in self.createMoleculeDict[molName][1]:
        self.localCreateMoleculeDict[molName][1].append((sequence,sequenceList[:]))
        
    sequences = self.localCreateMoleculeDict[self.molName][1]
    (self.origSequence,self.origSequenceList) = sequences[0]
          
  def apply(self):
    
    self.updateMolDict()

    return True
